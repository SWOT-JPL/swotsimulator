"""
@author <sylvain.herledan@oceandatalab.com>
@date 2018-12-26
"""

import os
import sys
import math
import time
import logging
import traceback
import multiprocessing

logger = logging.getLogger(__name__)
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.DEBUG)


class DyingOnError(Exception):
    """"""
    pass


class MultiprocessingError(Exception):
    """"""
    def __init__(self, exc):
        """"""
        self.exc = exc


class JobsManager():
    """"""
    def __init__(self, pool_size, status_updater, exc_fmt, err_fmt):
        """"""
        self._pool_size = pool_size
        self.manager = multiprocessing.Manager()
        self.msg_queue = self.manager.Queue()
        self.errors_queue = self.manager.Queue()
        self.pool = multiprocessing.get_context('spawn').Pool(pool_size)

        #self.format_exception = (exc_fmt,)
        self.format_exception = None
        if exc_fmt is not None:
            self.format_exception = staticmethod(exc_fmt)

        self.update_status = None
        if status_updater is not None:
            self.update_status = staticmethod(status_updater)

        self.format_error = None
        if err_fmt is not None:
            self.format_error = staticmethod(err_fmt)


    def show_errors(self):
        """"""
        if self.format_error is None:
            return

        while not self.errors_queue.empty():
            (pid, job_id, extra, exc) = self.errors_queue.get()
            error_str = self.format_error.__func__(pid, job_id, extra, exc)
            logger.error(error_str)
            logger.error('  {}'.format('  \n'.join(exc)))


    def _error_callback(self, exc):
        """"""
        logger.error(exc)
        raise MultiprocessingError(exc)



    def handle_message(self, status, msg, die_on_error, progress_bar):
        """die_on_error"""
        _ok = (msg[3] is None)
        if (progress_bar is True) and (self.update_status is not None):
            _ok = self.update_status.__func__(status, msg)
            #_ok = mod_tools.update_progress_multiproc(status, msg)
        if (_ok is False) and (die_on_error is True):
            # Kill all workers, show error and exit with status 1
            self.pool.terminate()
            self.show_errors()
            raise DyingOnError
        return _ok


    def submit_jobs(self, operation, jobs, die_on_error, progress_bar,
                    delay=0.5):
        """"""
        for j in jobs:
            j.append(self.errors_queue)
            j.append(self.msg_queue)
            j.append(self.format_exception.__func__)
            j.append(operation)

        # Distribute jobs between workers
        chunk_size = int(math.ceil(len(jobs) / self._pool_size))
        status = {}
        for n, w in enumerate(self.pool._pool):
            status[w.pid] = {'done': 0, 'total': 0, 'jobs': None, 'extra': ''}
            proc_jobs = jobs[n::self._pool_size]
            status[w.pid]['jobs'] = [j[0] for j in proc_jobs]
            status[w.pid]['total'] = len(proc_jobs)

        # Make some room for the progress bars
        if progress_bar is True:
            sys.stdout.write('\n' * self._pool_size)
            sys.stdout.flush()

        # Start jobs processing 
        tasks = self.pool.map_async(_operation_wrapper, jobs,
                                    chunksize=chunk_size,
                                    error_callback=self._error_callback)

        # Wait until all jobs have been executed
        ok = True
        while not tasks.ready():
            if not self.msg_queue.empty():
                msg = self.msg_queue.get()
                _ok = self.handle_message(status, msg, die_on_error,
                                          progress_bar)
                ok = ok and _ok
            time.sleep(delay)

        # Make sure all messages have been processed
        while not self.msg_queue.empty():
            msg = self.msg_queue.get()
            _ok = self.handle_message(status, msg, die_on_error, progress_bar)
            ok = ok and _ok

        # Flush stdout buffer to avoid partial output issues 
        sys.stdout.flush()

        # Wait for workers to release resources and synchronize with the main
        # process
        self.pool.close()
        self.pool.join()
        #"""
        return ok


def _operation_wrapper(*args, **kwargs):
    """"""
    _args = args[0]
    operation = _args.pop()
    format_exc = _args.pop()
    msg_queue = _args.pop()
    errors_queue = _args.pop()

    try:
        job_id = _args[0]
        operation(msg_queue, *_args, **kwargs)
    except:
        # Error sink
        exc = sys.exc_info()
        if format_exc is None:
            error_msg = traceback.format_exception(exc[0], exc[1], exc[2])
        else:
            error_msg = format_exc(exc)

        # Pass the error message to both the messages queue and the
        # errors queue
        msg_queue.put((os.getpid(), job_id, -1, error_msg))
        errors_queue.put((os.getpid(), job_id, -1, error_msg))
        return False

    return True
