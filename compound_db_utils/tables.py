import warnings

from sqlalchemy import Table, MetaData
from sqlalchemy import create_engine
from sqlalchemy.exc import SAWarning

from compound_db_utils import settings

_db = create_engine(
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
_metadata = MetaData(_db)

COMPOUNDS = None
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=SAWarning)
    COMPOUNDS = Table(
        settings.COMPOUNDS_TABLE
        , _metadata
        , autoload=True
        , include_columns=[
            'id'
            , 'unique_id'
            , 'smiles'
        ]
    )

DATASETS = Table(
    settings.DATASETS_TABLE
    , _metadata
    , autoload=True
)

BIOACTIVITIES = Table(
    settings.BIOACTIVITIES_TABLE
    , _metadata
    , autoload=True
)