import os
import sys
import functools
import subprocess
import platform
import textwrap

import pytest
import jaraco.envs
import path


IS_PYPY = '__pypy__' in sys.builtin_module_names


class VirtualEnv(jaraco.envs.VirtualEnv):
    name = '.env'
    # Some version of PyPy will import distutils on startup, implicitly
    # importing setuptools, and thus leading to BackendInvalid errors
    # when upgrading Setuptools. Bypass this behavior by avoiding the
    # early availability and need to upgrade.
    create_opts = ['--no-setuptools']

    def run(self, cmd, *args, **kwargs):
        cmd = [self.exe(cmd[0])] + cmd[1:]
        return subprocess.check_output(cmd, *args, cwd=self.root, **kwargs)


@pytest.fixture
def venv(tmp_path, tmp_src):
    env = VirtualEnv()
    env.root = path.Path(tmp_path / 'venv')
    env.req = str(tmp_src)
    return env.create()


def popen_text(call):
    """
    Augment the Popen call with the parameters to ensure unicode text.
    """
    return functools.partial(call, universal_newlines=True) \
        if sys.version_info < (3, 7) else functools.partial(call, text=True)


def win_sr(env):
    """
    On Windows, SYSTEMROOT must be present to avoid

    > Fatal Python error: _Py_HashRandomization_Init: failed to
    > get random numbers to initialize Python
    """
    if env is None:
        return
    if platform.system() == 'Windows':
        env['SYSTEMROOT'] = os.environ['SYSTEMROOT']
    return env


def find_distutils(venv, imports='distutils', env=None, **kwargs):
    py_cmd = 'import {imports}; print(distutils.__file__)'.format(**locals())
    cmd = ['python', '-c', py_cmd]
    return popen_text(venv.run)(cmd, env=win_sr(env), **kwargs)


def count_meta_path(venv, env=None):
    py_cmd = textwrap.dedent(
        """
        import sys
        is_distutils = lambda finder: finder.__class__.__name__ == "DistutilsMetaFinder"
        print(len(list(filter(is_distutils, sys.meta_path))))
        """)
    cmd = ['python', '-c', py_cmd]
    return int(popen_text(venv.run)(cmd, env=win_sr(env)))


def test_distutils_stdlib(venv):
    """
    Ensure stdlib distutils is used when appropriate.
    """
    env = dict(SETUPTOOLS_USE_DISTUTILS='stdlib')
    assert venv.name not in find_distutils(venv, env=env).split(os.sep)
    assert count_meta_path(venv, env=env) == 0


def test_distutils_local_with_setuptools(venv):
    """
    Ensure local distutils is used when appropriate.
    """
    env = dict(SETUPTOOLS_USE_DISTUTILS='local')
    loc = find_distutils(venv, imports='setuptools, distutils', env=env)
    assert venv.name in loc.split(os.sep)
    assert count_meta_path(venv, env=env) <= 1


@pytest.mark.xfail('IS_PYPY', reason='pypy imports distutils on startup')
def test_distutils_local(venv):
    """
    Even without importing, the setuptools-local copy of distutils is
    preferred.
    """
    env = dict(SETUPTOOLS_USE_DISTUTILS='local')
    assert venv.name in find_distutils(venv, env=env).split(os.sep)
    assert count_meta_path(venv, env=env) <= 1


def test_pip_import(venv):
    """
    Ensure pip can be imported.
    Regression test for #3002.
    """
    cmd = ['python', '-c', 'import pip']
    popen_text(venv.run)(cmd)
