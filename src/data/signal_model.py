import os
import pyabf
import neo
from loguru import logger
from src import settings

DEFAULT_STEP = (5, 0.1, 0.7)


class Signal:
    def __init__(self, filename, step=DEFAULT_STEP):
        self.filename = filename
        self.step = step
        self.report_folder = settings.REPORTS_DIR_RECORDING / filename.stem

        logger.info("Loading data...")
        self.file = self._read_file_by_neo(filename)
        self.create_report_folder(self.report_folder)
        self._read_file_by_pyabf(filename)

    @staticmethod
    def create_report_folder(report_folder):
        if os.path.exists(report_folder):
            pass
        else:
            os.mkdir(report_folder)

    def _read_file_by_neo(self, filename):
        logger.info("Reading file.")
        # return neo.io.AxonIO(filename)
        return neo.io.get_io(filename)

    def _read_file_by_pyabf(self, filename):
        self.abf = pyabf.ABF(filename)
        self.sweep_number = self.abf.sweepCount

        block = self.create_block()
        self.sampling_rate = block.segments[0].analogsignals[0].sampling_rate.item()
        self.raw_signal = block.segments[0].analogsignals[0].base.T
        self.times = [
            i / block.segments[0].analogsignals[0].sampling_rate.item()
            for i in range(len(self.raw_signal[0]))
        ]
        self.len_recording_time = len(self.raw_signal[0]) / self.sampling_rate
        self.sweep_length = int(self.abf.sweepIntervalSec)

    def create_block(self):
        return self.file.read_block()
