#"Cell 1 ctrl 4m", "Cell 2 ctrl 4m", "Cell 3 ctrl 4m", "Cell 4 ctrl 4m", "Cell 1 ctrl 8m", "Cell 2 ctrl 8m", "Cell 3 ctrl 8m", "Cell 4 ctrl 8m", "Cell 1 ad 4m", "Cell 2 ad 4m", "Cell 3 ad 4m", "Cell 4 ad 4m", "Cell 5 ad 4m", "Cell 6 ad 4m", "Cell 7 ad 4m", "Cell 1 ad 8m", "Cell 4 ad 8m", "Cell 6 ad 8m"

import uvicorn
import os
from loguru import logger

from src import settings
from fastapi import FastAPI, UploadFile, Form
from src.visualization.common_statistic import *

from src.data.signal_model import Signal
from src.features.signal_features import SignalFeatures
from src.tools.build_common_statistic_function import read_electric_cell_parameters
from src.visualization.report import make_report
import warnings
warnings.filterwarnings("ignore")

from pydantic import BaseModel
import typing

app = FastAPI()


class GroupAnalysisRequest(BaseModel):
    files: typing.List[str]
    report_folder: str
    colors: typing.Dict[str, typing.Tuple[int, ...]]

class RecordingAnalysisRequest(BaseModel):
    colors: typing.Dict[str, typing.Tuple[int, ...]]


STEP = (5, 0.1, 0.7)  # voltage, start, stop
COLORS = {
    "CTRL": (0, 30, "green"),
    "GABA1": (54, 72, "blue"),
    "GABA5": (100, 118, "navy"),
    "PTX": (162, 198, "red"),
}
# COLORS = {
#     "CTRL": [0, 30],
#     "GABA1": [54, 72],
#     "GABA5": [100, 118],
#     "PTX": [162, 198]
# }

def recording_analysis(filename, colors):
    signal = Signal(filename=filename, step=STEP)
    signal_features = SignalFeatures(signal=signal, colors=colors)
    signal_features.calculate_statistics()
    make_report(signal, signal_features)

    logger.info("Pipeline finished successfully.")


@app.post("/electrophysiology_recording_analysis/")
def electrophysiology_recording_analysis(file: UploadFile):
    if file is not None:
        contents = file.file.read()
        filename = (
            settings.RAW_DIR / file.filename
        )  # settings.RAW_DIR / "Cell 1 ctrl 4m.abf"
        with open(filename, "wb") as f:
            f.write(contents)
    else:
        filename = settings.RAW_DIR / "Cell 1 ctrl 4m.abf"

    #colors = {k: (v[0], v[1], settings.DEFAULT_COLORS[i]) for i, (k, v) in enumerate(colors.items())}
    recording_analysis(filename, COLORS)
    return {'status': 'OK'}

def build_common_statistic(report_folder, files, colors):
    electric_parameters_group_common_statistic(report_folder, files, colors)
    logger.info("Electric cell parameters common statistic is creating.")
    baseline_analysis_group_common_statistic(report_folder, files, colors)
    logger.info("Baseline parameters common statistic is creating.")
    event_parameters_group_common_statistic(report_folder, files, colors)
    logger.info("Event parameters common statistic is creating.")
    logger.info("Pipeline finished successfully.")

@app.post("/group_analysis/")
def group_analysis(group_request: GroupAnalysisRequest):
    files = group_request.files
    for filename in files:
        is_report_folder = settings.REPORTS_DIR_RECORDING / filename
        if not os.path.exists(is_report_folder):
            logger.warning(f"Report folder for file `{filename}` does not exist.\n" 
                           "Starting of single analysis pipeline.")
            recording_analysis(settings.RAW_DIR / f'{filename}.abf')

    report_folder = settings.REPORTS_DIR_GROUP / group_request.report_folder
    Signal.create_report_folder(report_folder)

    colors = {k: (v[0], v[1], settings.DEFAULT_COLORS[i]) for i, (k, v) in enumerate(group_request.colors.items())}
    build_common_statistic(report_folder, files, colors)


if __name__ == "__main__":
    # electrophysiology_recording_analysis(None)
    uvicorn.run("main:app", port=8000, log_level="info")
