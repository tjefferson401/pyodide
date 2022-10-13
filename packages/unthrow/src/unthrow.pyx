#cython: language_level=3
from cpython.object cimport PyObject,Py_SIZE
from cpython.ref cimport Py_INCREF,Py_XDECREF,Py_XINCREF,Py_CLEAR
from libc.string cimport memcpy
from cpython cimport array
import array
import collections
import sys,inspect,dis
# TJ IMPORTS
import time
# import json
# import ast
import copy
# END TJ IMPORTS

traceall=False

DEF PythonVersion=3.9

class Resumer:
    def __init__(self):
        self.finished=False
        self.resume_params=None
        self.resume_stack=None
        self.freq=0
        self.running=0
        # TJ EDIT BEGIN
        self.local_vars=None
        # self.first_line_run=True
        # TJ EDIT END

    def cancel(self):
        global __skip_stop
        self.resume_stack=None
        __skip_stop=False
        set_resume(0)

    def run_once(self,mainfn,args):
        global interrupts_enabled
        if self.resume_stack:
            _do_resume(<PyObject*>self.resume_stack) # --> _do_resume
            self.resume_stack=None
            # TJ EDIT BEGIN
            self.local_vars=None
            # TJ EDIT END
        self.running=1
        interrupts_enabled=1
        interrupt_call_level=-1
        # we call into the first level where interrupts should be
        interrupt_with_block_initial_level=PyEval_GetFrame().f_iblock+1
        self.finished=True

        set_interrupt_frequency(self.freq)

        try:
            mainfn(*args)
        except ResumableException as re:
            self.resume_params=re.parameter
            self.resume_stack=re.saved_frames
            # TJ EDIT BEGIN TODO
            self.local_vars=re.local_vars
            # TJ EDIT END
            re.with_traceback(None)
            self.finished=False
        finally:
            interrupts_enabled=0
        # stop interrupts
        global interrupts_enabled
        set_interrupt_frequency(0)
        self.running=0
        interrupts_enabled=0
        return self.finished

    def set_interrupt_frequency(self,freq):
        self.freq=freq
        if self.running:
            set_interrupt_frequency(self.freq)



# like a named tuple but mutable (so that we can zero things to avoid leaks)
class _SavedFrame:
    def __init__(self,*,locals_and_stack,lasti,code,block_stack,globals_if_different,slow_locals):
        dict=locals().copy() # put all the parameters to this into the dict
        del dict["self"]
        self._dictionary=dict

    def __getattr__(self,attr):
        if attr!="_dictionary" and attr in self._dictionary:
            return self._dictionary[attr]
        else:
            super().__getattribute__(attr)

    def __setattr__(self,attr,value):
        if attr!="_dictionary" and attr in self._dictionary:
            self._dictionary[attr]=value
        else:
            super().__setattr__(attr,value)

    def __repr__(self):
        return str(self._dictionary)

class _PythonNULL(object):
    pass

def test_iter(maxVal):
    for c in range(maxVal):
        print("TI:",c)
        yield c


cdef int  __skip_stop=False

cdef extern from "frameobject.h":
    ctypedef struct PyTryBlock:
        int b_type
        int b_handler
        int b_level


    cdef struct _frame:
        _frame * f_back
        PyCodeObject* f_code
        int f_lasti
        # TJ ADD
        int f_lineno
        # TJ ADD END
        char f_executing
        PyObject *f_locals
        PyObject *f_globals
        PyObject *f_builtins

        PyObject **f_stacktop
        PyObject **f_valuestack
        PyObject **f_localsplus
        PyTryBlock f_blockstack[1] # this is actually sized by a constant, but cython
                                   # doesn't redeclare it so we can just put 1 in
        int f_iblock
        int f_trace_opcodes;


    ctypedef _frame PyFrameObject

    cdef PyFrameObject* PyEval_GetFrame()
    cdef void PyFrame_FastToLocals(PyFrameObject* frame)
    cdef void PyFrame_LocalsToFast(PyFrameObject* frame,int clear)
    cdef PyObject* PyObject_Call(PyObject*obj,PyObject*args,PyObject*kwArgs)


cdef extern from "opcode.h":
    IF PythonVersion==3.8:
        cdef enum _opcodes:
            SETUP_FINALLY
            BEGIN_FINALLY
    ELSE:
        # 3.9 or newer
        cdef enum _opcodes:
            SETUP_FINALLY
            DUP_TOP
            POP_TOP
            SETUP_WITH
            SETUP_ASYNC_WITH


cdef extern from "Python.h":
    ctypedef struct PyCodeObject:
        PyObject* co_code
        PyObject* co_name
        PyObject* co_filename
        int co_flags

    cdef PyObject* Py_True
    ctypedef int (*Py_tracefunc)(PyObject *obj, PyFrameObject *frame, int what, PyObject *arg)
    cdef enum _trace_types:
        PyTrace_LINE
        PyTrace_CALL
        PyTrace_C_CALL
        PyTrace_OPCODE
        PyTrace_RETURN

    cdef int PyCompile_OpcodeStackEffectWithJump(int,int,int)
    cdef char* PyBytes_AsString(PyObject *o)
    cdef PyObject* PyBytes_FromStringAndSize(const char *v, Py_ssize_t len)
    cdef Py_ssize_t PyObject_Length(PyObject *o)
    cdef PyObject* PyDict_GetItem(PyObject*,PyObject*)
    cdef int PyDict_SetItemString(PyObject*dict,const char*key,PyObject* val)
    cdef int PyDict_SetItem(PyObject*dict,PyObject*key,PyObject* val)
    cdef int PyDict_Next(PyObject *p, Py_ssize_t *ppos, PyObject **pkey, PyObject **pvalue)

    cdef int PyCode_Check(PyObject*)

    cdef int PyEval_SetTrace(Py_tracefunc fn,PyObject*self)
    cdef int Py_AddPendingCall(int (*fn)(void*),void*arg)
    cdef void PyErr_SetObject(PyObject*type,PyObject*val)
    cdef PyObject* PyUnicode_FromString(char*str);
    cdef PyObject* PyExc_ValueError

    cdef PyObject* PyDict_New()

    cdef void PyErr_SetString(PyObject*,char *)
    cdef PyObject* PyErr_NewException(const char *name, PyObject *base, PyObject *dict)
    cdef PyObject* PyErr_Occurred()
    cdef PyObject* PyTuple_New(Py_ssize_t size)

    cdef int Py_REFCNT(PyObject* ob)

    cdef PyObject* PyList_GetItem(PyObject*,Py_ssize_t pos)

    cdef PyObject* PyObject_Type(PyObject*ob)
    cdef int PyGen_Check(PyObject* ob)
    ctypedef struct PyGenObject:
        PyFrameObject *gi_frame

cdef PyObject* ResumableExceptionClass

ResumableExceptionClass=<PyObject*>PyErr_NewException("unthrow.ResumableException",NULL,NULL)
ResumableException=<object>ResumableExceptionClass

cdef _get_stack_pos(object code,int target,int before):
    if target<0: # start of fn, empty stack
        return 0
    cdef int no_jump
    cdef int yes_jump
    cdef int opcode
    cdef int arg
    stack_levels={}
    jump_levels={}
    cur_stack=0
    for i in dis.get_instructions(code):
        offset=i.offset
        argval=i.argval
        if i.arg:
            arg=int(i.arg)
        else:
            arg=0
        opcode=int(i.opcode)
        if offset in jump_levels:
            cur_stack=jump_levels[offset]
        no_jump=PyCompile_OpcodeStackEffectWithJump(opcode,arg,0)
        if opcode in dis.hasjabs or opcode in dis.hasjrel:
            # a jump - mark the stack level at jump target
            yes_jump=PyCompile_OpcodeStackEffectWithJump(opcode,arg,1)
            if not argval in jump_levels:
                jump_levels[argval]=cur_stack+yes_jump
        if before==1:
            stack_levels[offset]=cur_stack
        cur_stack+=no_jump
        if before==0:
            stack_levels[offset]=cur_stack
#        print(target,offset,i.opname,argval,cur_stack)
#    print(stack_levels[target])
    return stack_levels[target]

cdef get_stack_pos_before(object code,int target):
    return _get_stack_pos(code,target,1)

cdef get_stack_pos_after(object code,int target):
    return _get_stack_pos(code,target,0)

cdef get_line_start(object code,int target):
    lineStart=target
    for i in dis.get_instructions(code):
        if i.offset>=target:
            break
        if i.starts_line:
            lineStart=i.offset
#    print("Line start for ",target,"=",lineStart,end="")
#    dis.dis(code)
    return lineStart

# TJ BEGIN
cdef object step_info(PyFrameObject* source_frame):
    slow_locals=None
    lineno = 0
    if source_frame.f_code:
        slow_locals=(<object>source_frame).f_locals.copy()
        lineno = (<object>source_frame).f_lineno
        code_name = (<object>source_frame).f_code.co_name
    return (slow_locals, lineno, (<object>source_frame).f_code.co_name)

# TJ END

cdef object save_frame(PyFrameObject* source_frame,from_interrupt):
    print("save_frame")
    cdef PyObject *localPtr;
    if (source_frame.f_code.co_flags & inspect.CO_OPTIMIZED)!=0:
        PyFrame_LocalsToFast(source_frame,0)
    blockstack=[]
    # last instruction called
    lasti=source_frame.f_lasti
    # get our value stack size from the code instructions
    if from_interrupt:
        # in an interrupt top level frame, the function is at the start of an instruction
        # and the next instruction is about to run
        lasti-=2
        if source_frame.f_stacktop!=NULL:
            our_stacksize=source_frame.f_stacktop - source_frame.f_valuestack
        else:
            our_stacksize=get_stack_pos_before(<object>source_frame.f_code,lasti+2)
    else:
        # anywhere else, the frame is half way through a call and lasti points
        # at just the call. But we still have the call args pushed on stack,
        # lasti = lasti-2 to make the next thing be the call again
        # i.e. we want to run this instruction twice
        lasti=get_line_start(<object>source_frame.f_code,lasti)
        lasti-=2
        our_stacksize=get_stack_pos_before(<object>source_frame.f_code,lasti+2)
    our_localsize=<int>(source_frame.f_valuestack-source_frame.f_localsplus);
    code_obj=(<object>source_frame).f_code.co_code
    # grab everything off the locals and value stack
    valuestack=[]
    for c in range(our_stacksize+our_localsize):
        if <int>source_frame.f_localsplus[c] ==0:
            valuestack.append(_PythonNULL())
        else:
            localObj=<object>source_frame.f_localsplus[c]
            valuestack.append(localObj)
    blockstack=[]
    for c in range(source_frame.f_iblock):
        blockstack.append(source_frame.f_blockstack[c])
    globals_if_different=None
    # if we are in a module (or a weird exec call)
    # then save module level globals also
    # n.b. this is a shallow copy, so it probably will only work nicely
    # within the same interpreter - I don't think it is going to be
    # persistible
    if source_frame.f_back!=NULL and source_frame.f_back.f_globals != source_frame.f_globals:
        # ORIGINAL CODE BELOW:
        # globals_if_different=(<object>source_frame).f_globals.copy()
        # print(type((<object>source_frame).f_globals))
        # TJ ADD CODE ADDED Thu Sep 22 10:37 PM
        # sys.setrecursionlimit(5000)
        print("save_frame: Copying globals if different iterively")
        globals_if_different=iter_deepcopy((<object>source_frame).f_globals)
        # sys.setrecursionlimit(1000)
        # THIS DID NOT WORK. "Error: Internal error: Argument 'undefined' to hiwire.get_value is falsy (but error indicator is not set)."
        # globals_if_different= ast.parse(ast.dump((<object>source_frame).f_globals))
        # SECOND ATTEMPT
        # deepcopy((<object>source_frame).f_globals)
        # TJ END CODE
        # don't save the builtins dict
        if "__builtins__" in globals_if_different:
            del globals_if_different["__builtins__"]

    # if we are in a non-optimized frame, i.e. without fast locals, we need to copy the locals dict
    slow_locals=None
    if (source_frame.f_code.co_flags & inspect.CO_OPTIMIZED)==0 and source_frame.f_locals!=source_frame.f_globals:
        print("save_frame: Copying slow locals")
        slow_locals=(<object>source_frame).f_locals.copy()


    return _SavedFrame(locals_and_stack=valuestack,code=code_obj,lasti=lasti,block_stack=blockstack,globals_if_different=globals_if_different,slow_locals=slow_locals)

cdef void restore_saved_frame(PyFrameObject* target_frame,saved_frame: _SavedFrame):
    print("restore_saved_frame")
    cdef PyObject* tmpObject
    cdef PyObject* borrowed_list_item;
    if (target_frame.f_code.co_flags & inspect.CO_OPTIMIZED)!=0:
        PyFrame_LocalsToFast(target_frame,1)
    # last instruction
    target_frame.f_lasti=saved_frame.lasti
    # check code is the same
    if (<object>target_frame).f_code.co_code!=saved_frame.code:
        print("Trying to restore wrong frame")
        return
    # restore locals and stack
    # n.b. we currently know we're at the start of the stack, but we do need
    # to free any locals that are currently non-null
    num_locals=<int>(target_frame.f_valuestack-target_frame.f_localsplus)
    num_oldstack=<int>(target_frame.f_stacktop-target_frame.f_valuestack)
    for c in range(len(saved_frame.locals_and_stack)):
        print("Looping through locals and stack")
        if c<num_locals:
            tmpObject=target_frame.f_localsplus[c]
        else:
            tmpObject=NULL
        borrowed_list_item=PyList_GetItem(<PyObject*>saved_frame.locals_and_stack,c)
        if PyObject_Type(borrowed_list_item)==<PyObject*>_PythonNULL:
            target_frame.f_localsplus[c]=NULL
        else:
            target_frame.f_localsplus[c]=borrowed_list_item
            Py_XINCREF(target_frame.f_localsplus[c])
        if tmpObject!=NULL:
            Py_XDECREF(tmpObject)



    target_frame.f_stacktop=&target_frame.f_localsplus[len(saved_frame.locals_and_stack)]
    saved_frame.locals_and_stack=[]
    # restore block stack
    for c,x in enumerate(saved_frame.block_stack):
        target_frame.f_blockstack[c]=x
    target_frame.f_iblock=len(saved_frame.block_stack)
    saved_frame.block_stack=[]
    # if globals are different (i.e. in a module), then restore them

    cdef PyObject *key
    cdef PyObject *srcValue
    cdef PyObject *targetValue
    cdef Py_ssize_t pos = 0

    if saved_frame.globals_if_different!=None:
        while PyDict_Next(<PyObject*>(saved_frame.globals_if_different), &pos, &key, &srcValue)!=0:
            targetValue=PyDict_GetItem(target_frame.f_globals,key)
            # add any new keys
            set_it=False
            if targetValue==NULL:
                set_it=True
            elif targetValue!=srcValue:
                # if they point to a different object set it
                set_it=True
            if set_it:
                #print("RESTORING:",<object>key,<object>srcValue)
                PyDict_SetItem(target_frame.f_globals,key,srcValue)
    saved_frame.globals_if_different=None

    if (target_frame.f_code.co_flags & inspect.CO_OPTIMIZED)==0 and saved_frame.slow_locals!=None:
        while PyDict_Next(<PyObject*>(saved_frame.slow_locals), &pos, &key, &srcValue)!=0:
            targetValue=PyDict_GetItem(target_frame.f_locals,key)
            # add any new keys
            set_it=False
            if targetValue==NULL:
                set_it=True
            elif targetValue!=srcValue:
                # if they point to a different object set it
                set_it=True
            if set_it:
                #print("RESTORING:",<object>key,<object>srcValue)
                PyDict_SetItem(target_frame.f_locals,key,srcValue)
    saved_frame.slow_locals=None
    del saved_frame
    print("restore_saved_frame ENDING")


# the parent frame where interrupts were started.
# if this is not in the call stack then we don't throw an exception
cdef int interrupts_enabled=1
cdef int interrupt_counter=0
cdef int interrupt_frequency=0 # every N instructions call out to JS
cdef int in_resume=0
cdef int interrupt_with_block_initial_level=0
cdef int interrupt_call_level=0
cdef int interrupt_with_level=-1

cdef void set_resume(int is_resuming):
    global in_resume,resume_state,interrupt_frequency,_trace_obj
    resume_state=0
    in_resume=is_resuming
    if in_resume==0 and interrupt_frequency==0:
        print("SETTING TRACE IF")
        PyEval_SetTrace(NULL,NULL)
        print("AFTER TRACE if")
    else:
        print("SETTING TRACE ELSE")
        PyEval_SetTrace(&_c_trace_fn,NULL)
        print("AFTER TRACE ELSE")

cdef void set_interrupt_frequency(int freq):
    global in_resume,resume_state,interrupt_frequency,interrupt_counter,_trace_obj
    if interrupt_frequency!=freq:
        interrupt_counter=0
        interrupt_frequency=freq
        if in_resume==0 and interrupt_frequency==0:
            PyEval_SetTrace(NULL,NULL)
        else:
            PyEval_SetTrace(&_c_trace_fn,NULL)

# check if this frame is currently within a with or finally block
cdef int _check_blocks(PyFrameObject* frame):
    frameCode=PyBytes_AsString(frame.f_code.co_code)
    for c in range(frame.f_iblock):
        # inside a with block
        if frame.f_blockstack[c].b_type==SETUP_WITH or frame.f_blockstack[c].b_type==SETUP_ASYNC_WITH:
            return True
        # inside a block which is a try, finally block
        if (frame.f_blockstack[c].b_type)==SETUP_FINALLY: # with block is basically try, finally
            # SETUP_FINALLY also sets up exception catches, an actual finally block
            # on 3.8
            # is a SETUP_FINALLY that points to POP_BLOCK, BEGIN_FINALLY
            # on 3.9 it is a setupfinally that points to
            # something other than DUP_TOP or POP_TOP
            # we only care about exception blocks that would be called when our exception
            # comes out. nb. catch-all exception blocks will also get called
# 3.8 code
#                                if frameCode[frame.f_blockstack[c].b_handler-2]==BEGIN_FINALLY:
#                                   # there is a finally block
#                                    interrupt_with_level=interrupt_call_level
            if frameCode[frame.f_blockstack[c].b_handler]!=DUP_TOP and frameCode[frame.f_blockstack[c].b_handler]!=POP_TOP:
               # there is a finally block
                return True
            else:
                return False

cdef int _c_trace_fn(PyObject *self, PyFrameObject *frame,
                 int what, PyObject *arg):
    global interrupt_frequency,interrupt_counter,interrupts_enabled,interrupt_call_level,interrupt_with_level
    # ADD TJ
    # ADD TJ END
    if in_resume:
        # in resume call, ignore interrupts
        if what==PyTrace_CALL:
            _resume_frame(_resume_list,frame)
    elif interrupts_enabled==1:
        if what==PyTrace_CALL:
            # check if this call is enter or exit of a with
            if <object>(frame.f_code.co_name)=="__enter__" or <object>(frame.f_code.co_name)=="__exit__":
                if interrupt_with_level==-1:
                    interrupt_with_level=interrupt_call_level
            # and don't interrupt in import machinery
            elif (<object>(frame.f_code.co_filename)).find("importlib.")!=-1:
                if interrupt_with_level==-1:
                    interrupt_with_level=interrupt_call_level
            # and check if the parent is within a finally block etc.
            elif interrupt_call_level!=-1 and _check_blocks(frame.f_back):
                if interrupt_with_level==-1:
                    interrupt_with_level=interrupt_call_level
            interrupt_call_level+=1
        if what==PyTrace_RETURN:
            if interrupt_call_level<interrupt_with_level:
                interrupt_with_level=-1
            interrupt_call_level-=1
            if interrupt_call_level<0:
                # outside the scope where we were called, stop tracing now
                set_interrupt_frequency(0)
        if interrupt_frequency!=0 and what==PyTrace_LINE:
            interrupt_counter+=1
            # only throw interrupt if we are not inside a with block,
            # as we can't restore context managers (files etc)
            if interrupt_counter>=interrupt_frequency:
                # check that this frames isn't within a with or finally block
                # (parent frames are checked on call)
                if interrupt_with_level==-1 or interrupt_with_level>=interrupt_call_level:
                    interrupt_with_level=-1
                    if _check_blocks(frame):
                      interrupt_with_level=interrupt_call_level
                interrupts_enabled=0
                if interrupt_with_level==-1:
                    # throw interrupt exception
                    interrupt_counter=0
                    make_interrupt(<void*>self,frame)
                    return 1 # need to return 1 to signal error or else our exception
                             # gets cleaned up by cpython
                else:
                    interrupts_enabled=1

    return 0

cdef make_interrupt(void* arg,PyFrameObject*frame):
    interrupts_enabled=0
    cdef PyObject *rex
    cdef PyObject* rex_type
    cdef char*str="interrupt"
    cdef PyObject* msg=PyDict_New()
    PyDict_SetItemString(msg,"interrupt",Py_True)
    try:
        rex=make_resumable_exception(msg,frame)
        rex_type=<PyObject*>rex.ob_type
        PyErr_SetObject(<PyObject*>rex_type,rex)
        Py_XDECREF(rex)
    except Exception as err:
        print("PROBLEM IN CONSTRUCT REX",err)
    return 0


cdef PyObject* make_resumable_exception(PyObject* msg,PyFrameObject* frame):
    cdef PyObject* exc=PyObject_Call(ResumableExceptionClass,PyTuple_New(0),NULL)
    (<object>exc).parameter=<object>msg;
    (<object>exc).saved_frames=[]
    # TJ ADD BEGIN
    (<object>exc).local_vars=[]
    #_save_stack((<object>exc).saved_frames,frame)
    _save_stack((<object>exc).saved_frames,frame,(<object>exc).local_vars)
    # TJ ADD END
    return exc



# store locals and stack for all the things, while we're still in a call, before exception is thrown
# and calling stack objects are dereferenced etc.
    # TJ ADD BEGIN
cdef _save_stack(object saved_frames,PyFrameObject* cFrame, object local_vars):
    from_interrupt=False
    if cFrame==NULL:
        cFrame=PyEval_GetFrame()
    else:
        from_interrupt=True
    while cFrame!=NULL:
        # TODO: SAVE local_vars
        # TJ BEGIN
        local_vars.append(step_info(cFrame))
        # Tj END
        saved_frames.append(save_frame(cFrame,from_interrupt=from_interrupt))
        from_interrupt=False # only the top frame of an interrupt is different
        cFrame=cFrame.f_back

cdef void _do_resume(PyObject* c_saved_frames):
    # c_saved_frames = resumer.resume_stack
    # pulls in globals skip_stop, resume_list
    global __skip_stop,_resume_list
    # Increments resume_stack if it exists
    Py_XINCREF(c_saved_frames)
    # Sets saved_frames to copy of c_...
    saved_frames=<object>c_saved_frames
    # creates new object cFrame
    cdef PyFrameObject* cFrame;
    # Sets cFrame to the current thread state’s frame
    cFrame=PyEval_GetFrame()
    # init allStack
    allStack=[]
    # while cFrame is not null, add cFrame back to allStack
    while cFrame!=NULL:
        allStack.append(<object>cFrame)
        cFrame=cFrame.f_back
    # loop through allStack from back to front
    for frame in reversed(allStack):
        # if laset el of saved_frame.code = this frame, pop this frame
        if saved_frames[-1].code==frame.f_code.co_code:
            saved_frames.pop()
    # reset allStack and increment Saved frames
    allStack=[]
    Py_XINCREF(c_saved_frames)
    # Set resume list
    _resume_list=c_saved_frames
    set_resume(1)

cdef void _resume_frame(PyObject* c_saved_frames,PyFrameObject* c_frame):
    global interrupts_enabled
    frame=<object>c_frame
    saved_frames=<object>c_saved_frames
    if len(saved_frames)>0:
        resumeFrame=saved_frames[-1]
        if frame.f_code.co_code==resumeFrame.code:
            saved_frames.pop()
#            print("MATCH FRAME:",frame,resumeFrame)
            restore_saved_frame(c_frame,resumeFrame)
            if len(saved_frames)==0:
                print("Resuming")
                set_resume(0)
                print("AFTER SET RESUME")
                _resume_list=NULL;
                interrupts_enabled=1

# this is a c function so it doesn't get traced into
def stop(msg):
    global __skip_stop,interrupts_enabled
    if __skip_stop:
        #print("RESUMING - enable interrupts")
        __skip_stop=False
        interrupts_enabled=1
    else:
        __skip_stop=True
        # disable interrupts until we return from this
        interrupts_enabled=0
        #print("STOPPING NOW - disabled interrupts",<object>msg)
        rex=make_resumable_exception(<PyObject*>msg,NULL)
        objRex=<object>rex
        Py_XDECREF(rex)
        raise objRex

cdef PyObject* _resume_list=NULL





# run_once
#   _do_resume



def iter_deepcopy(obj):
    root = dict()
    stack = [(obj, root)]
    while stack:
        (original, deep), *stack = stack
        if isinstance(original, dict):
            for key in original:
                i = original[key]
                if isinstance(i, list):
                    p = (i, [])
                    deep[key] = p[1]
                    stack.append(p)
                elif isinstance(i, dict):
                    p = (i, dict())
                    deep[key] = p[1]
                    stack.append(p)
                else:
                    deep[key] = i
        elif isinstance(original, list):
            for i in original:
                if isinstance(i, list):
                    p = (i, [])
                    deep.append(p[1])
                    stack.append(p)
                elif isinstance(i, dict):
                    p = (i, dict())
                    deep.append(p[1])
                    stack.append(p)
                else:
                    deep.append(i)

        else:
            return;

    return root
