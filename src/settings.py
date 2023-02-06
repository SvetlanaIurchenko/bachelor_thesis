import pathlib
import typing

# project dirs
PROJECT_DIR: typing.Final[pathlib.Path] = pathlib.Path(__file__).parent.parent
DATA_DIR: typing.Final[pathlib.Path] = PROJECT_DIR / "data"
REPORTS_DIR: typing.Final[pathlib.Path] = PROJECT_DIR / "reports"
REPORTS_DIR_RECORDING: typing.Final[pathlib.Path] = REPORTS_DIR / "recordings"
REPORTS_DIR_GROUP: typing.Final[pathlib.Path] = REPORTS_DIR / "groups"
RAW_DIR: typing.Final[pathlib.Path] = DATA_DIR / "raw"
PROCESSED_DIR: typing.Final[pathlib.Path] = DATA_DIR / "processed"

# project paths
TEST_FILEPATH: typing.Final[pathlib.Path] = RAW_DIR / "test.abf"


DEFAULT_COLORS: typing.List[str] = ["green", "blue", "navy", "red"]
