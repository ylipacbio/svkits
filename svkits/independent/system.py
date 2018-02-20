"""Utilities which interact with the O/S
"""
import commands
import contextlib
import time
import datetime
import logging
import os
import os.path as op
import pipes
import subprocess

LOG = logging.getLogger()

@contextlib.contextmanager
def cd(newdir):
    # https://stackoverflow.com/questions/431684/how-do-i-cd-in-python
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def execute(cmd, errmsg="", errcls=RuntimeError):
    """Execute command and check exit code. If exit code is not
    0, raise RuntimeError with errmsg.
    """
    LOG.info('execute({!r})'.format(cmd))
    msg = ''
    code, out = commands.getstatusoutput(cmd)
    if code:
        msg = '{} <- `{}`\n{}'.format(code, cmd, out)
        raise errcls(msg)
    return out

def execute_as_bash(f):
    execute('chmod +x {f} && bash {f}'.format(f=f))

def quoted_abspath(f):
    """Return absolute path of f, shell-quoted if necessary.
    Expand ~ if any.
    """
    return pipes.quote(os.path.abspath(os.path.expanduser(f)))


def realpath(f):
    """Return absolute, user expanded path."""
    # copied from pbtranscript.Utils, like much here
    if f is None: # remove this may fail `pbtranscript` and `ice_*` scripts
        return None
    return op.abspath(op.expanduser(f))


def mkdir(d):
    if not os.path.isdir(d):
        LOG.info('mkdir -p {}'.format(d))
        os.makedirs(d)

def lnabs(src, dst):
    """if src and dst are identical, pass. Otherwise, create dst, a soft
    symbolic link pointing to abspath(src)."""
    if os.path.realpath(src) != os.path.realpath(dst):
        if op.exists(dst) or op.lexists(dst):
            os.remove(dst)
        logging.debug("Creating a symbolic link {dst} pointing to {src}".
                      format(dst=pipes.quote(dst), src=pipes.quote(src)))
        os.symlink(os.path.abspath(os.path.expanduser(src)), dst)

def mv(src, dst):
    """move src file to dst"""
    if os.path.realpath(src) != os.path.realpath(dst):
        execute('mv {} {}'.format(pipes.quote(src), pipes.quote(dst)))

def rmpath(path):
    """Remove a file or a directory"""
    if op.exists(path):
        execute("chmod -R +w {p} && rm -rf {p}".format(p=pipes.quote(path)))

def rmpath_without_chmod(path):
    """Remove a file or a directory assuming having write access to all files"""
    if op.exists(path):
        execute('rm -rf {p}'.format(p=pipes.quote(path)))

def rm(path):
    # TODO: merge rmpath and rm
    rmpath(path)

def touch(path):
    """touch a file."""
    if op.exists(path):
        os.utime(path, None)
    else:
        open(path, 'a').close()



def mknewdir(path):
    """Create a new directory if it does not pre-exist,
    otherwise, delete it and then re-create it."""
    rmpath(path)
    #shutil.rmtree(path)
    os.makedirs(path)


# Copied from pbcore.utils.Process
def backticks(cmd, merge_stderr=True):
    """
    Simulates the perl backticks (``) command with error-handling support
    Returns ( command output as sequence of strings, error code, error message )
    """
    LOG.info('backticks({!r}, {})'.format(
        cmd, merge_stderr))
    if merge_stderr:
        _stderr = subprocess.STDOUT
    else:
        _stderr = subprocess.PIPE

    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=_stderr,
                         close_fds=True)

    out = [l[:-1] for l in p.stdout.readlines()]

    p.stdout.close()
    if not merge_stderr:
        p.stderr.close()

    # need to allow process to terminate
    p.wait()

    errCode = p.returncode and p.returncode or 0
    if p.returncode > 0:
        errorMessage = os.linesep.join(out)
        output = []
    else:
        errorMessage = ''
        output = out

    return output, errCode, errorMessage


def nfs_exists(fn):
    """Detect whether a NFS file or a directory exists or not.
    In rare cases, a file f, which is created by a node X,
    may not be detected by node Y immediately due to NFS latency.
    In even more rare cases, node Y's shell script may be able
    to detect f's existence while Y's python can not, because of
    a python file system cache latency. So in order to eliminate
    this problem, first call 'ls' to trigger mount, then detect
    existence of fn using python.open(fn, 'r'), and try twice
    before finally give up.

    This script should return True, if fn is either
    an existing file or an existing directory created before
    nfs_exists is called. However, this script may return an
    arbitrary value, if fn is created or deleted at the same
    time or after nfs_exists is called.
    """
    if op.exists(fn):
        return True
    # Call ls just to trigger mount, don't trust the return value.
    _o, _c, _m = backticks("ls {f}".format(f=fn))

    ERROR_NO_SUCH_FILE_OR_DIRECTORY = 2
    ERROR_IS_DIRECTORY = 21
    # Try to open fn with read-only mode
    def check_once(fn):
        """Check file existance once."""
        try:
            with open(fn, 'r') as _reader:
                pass
            return True # fn is an existing file
        except IOError as err:
            if err.errno == ERROR_NO_SUCH_FILE_OR_DIRECTORY:
                return False # fn does not exist
            elif err.errno == ERROR_IS_DIRECTORY:
                return True # fn is an existing directory
            else:
                return False # other IOErrors

    # Retry every 10 seconds and max at 60 seoncds
    retry = 0
    OK = check_once(fn)
    while not OK and retry < 6:
        LOG.debug('Waiting for NFS file or directory {!r}'.format(fn))
        time.sleep(10)
        retry += 1
        OK = check_once(fn)
    return OK


def assert_nfs_exists(i_fn, desc):
    """Assert input file or dir must exist, """
    assert nfs_exists(i_fn), "{d} {f!r} must exists!".format(d=desc, f=i_fn)

def actually_now():
    """Return string of current time.
    """
    return datetime.datetime.now()

def tywmc():
    """the year we make contact
    '2001-01-01 00:00:00'
    """
    return datetime.datetime(2001, 1, 1, 0, 0, 0, 0)

def now_str():
    return str(now()).split(".")[0]

now = actually_now # can be replace by --log-test
