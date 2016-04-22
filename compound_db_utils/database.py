import warnings

from sqlalchemy import Table, MetaData
from sqlalchemy import create_engine
from sqlalchemy.exc import SAWarning

from compound_db_utils import settings

def fetch_all():
    DB_CONNECTION = create_engine(
        'postgresql+{0}://{1}:{2}@{3}:{4}/{5}'.format(
            settings.DATABASE['ENGINE']
            , settings.DATABASE['USER']
            , settings.DATABASE['PASSWORD']
            , settings.DATABASE['HOST']
            , settings.DATABASE['PORT']
            , settings.DATABASE['NAME']
        )
        , echo = False # Try changing this to True and see what happens
    )
    _metadata = MetaData(DB_CONNECTION)

    TB_COMPOUNDS = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=SAWarning)
        TB_COMPOUNDS = Table(
            settings.COMPOUNDS_TABLE
            , _metadata
            , autoload=True
            , include_columns=settings.COMPOUNDS_TABLE_COLS
        )

    TB_DATASETS = Table(
        settings.DATASETS_TABLE
        , _metadata
        , autoload=True
        , include_columns=settings.DATASETS_TABLE_COLS
    )

    TB_BIOACTIVITIES = Table(
        settings.BIOACTIVITIES_TABLE
        , _metadata
        , autoload=True
        , include_columns=settings.BIOACTIVITIES_TABLE_COLS
    )

    return DB_CONNECTION, TB_COMPOUNDS, TB_DATASETS, TB_BIOACTIVITIES