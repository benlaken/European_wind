import numpy as np
import pandas as pd
import datetime as dt
from calendar import monthrange
from dateutil.relativedelta import relativedelta

# Read the Sato index from the original file, convert the time index to
# datetime objects, and turn it all into a pandas dataframe. Pickle it 
# for convienent read in later.

def read_textfile(data_path):
    """ Read data from file with a generator"""
    with open(data_path) as f:
        for line in f:
            yield line

def decimal2date(atime):
    """
    Convert fractional time to pd.datetime, with the day
    as the last day of the month.
    atime is the decimal date (i.e. year.fraction) float.
    """
    assert  type(atime) == float
    year = int(atime)
    remainder = atime - year    
    boy = dt.datetime(year, 1, 1)
    eoy = dt.datetime(year + 1, 1, 1)
    seconds = remainder * (eoy - boy).total_seconds()
    mm = (boy + dt.timedelta(seconds=seconds)).month
    dd = monthrange(year, mm)[1]
    return pd.datetime(year,mm,dd).date()
                
sato_file = list(read_textfile("Data/Sato.txt"))
#print("Extracting Sato data for {0:s}".format(sato_file[3]))

sato_dates = []
nhemi_aod = []
for n, line in enumerate(sato_file[4:]):
    tmp_date, b, nhemi , d = [float(entry) for entry in line.split()]
    sato_dates.append(tmp_date)
    nhemi_aod.append(nhemi)

sato_dtindx = [decimal2date(sato_date) for sato_date in sato_dates] 

sato_df = pd.DataFrame(nhemi_aod,index=sato_dtindx,columns=['NHemi_AOD'])

sato_df.to_pickle("Data/sato_index")