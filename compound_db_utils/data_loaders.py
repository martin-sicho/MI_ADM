import sys

from compound_db_utils.tables import COMPOUNDS, DATASETS, BIOACTIVITIES

def main(args):
    activities = COMPOUNDS.select()
    for row in activities.execute():
        print(row.unique_id)



if __name__ == "__main__":
    exit(main(sys.argv))


