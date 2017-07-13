from sys import stderr, stdout
from logging import getLogger, DEBUG, Logger, Formatter, StreamHandler, ERROR, WARNING, INFO
from pprint import PrettyPrinter

LEVEL_FOR_VERBOSITY = dict(enumerate([ERROR, WARNING, INFO, DEBUG] + [DEBUG] * 10))

DEFAULT_FORMATTER = Formatter('%(asctime)s -[%(levelname)s]: %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')

USE_VERBOSITY_FILTER = True

class Verbosity_Filter(object):
    def __init__(self, verbosity):
        self.level = LEVEL_FOR_VERBOSITY[verbosity]

    def filter(self, log_record):
        return log_record.levelno >= self.level

def get_log(unique_id: str, verbosity: int = 0) -> Logger:
    log = getLogger(unique_id)
    if len(log.handlers) == 0:
        #fh = FileHandler('/tmp/log')
        #fh.setLevel(DEBUG)
        ch = StreamHandler(stream=stderr)
        ch.setLevel(DEBUG)
        #fh.setFormatter(DEFAULT_FORMATTER)
        ch.setFormatter(DEFAULT_FORMATTER)
        #log.addHandler(fh)
        log.addHandler(ch)
        log.setLevel(DEBUG)

        if USE_VERBOSITY_FILTER:
            log.addFilter(Verbosity_Filter(verbosity))

    else:
        pass

    return log

log = get_log('1', DEBUG)

pp = PrettyPrinter(indent=2)
pformat = pp.pformat
