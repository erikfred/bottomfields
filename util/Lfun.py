"""

Modified from Parker MacCready, School of Oceanography, University of Washington
https://github.com/parkermac

Functions for establishing file system, and other mechanics

"""
import os
from datetime import datetime, timedelta
import csv
import subprocess
import shutil

# this bit of magic lets us know where this program lives
# and so allows us to find get_lo_info.sh and lo_info.csv
alp = os.path.dirname(os.path.realpath(__file__))
alp_user = alp.replace('LiveOcean','LiveOcean_user')
upath = alp_user.replace('/alpha','')

def Lstart(gridname='BLANK', tag='BLANK'):
    """
    This is to set environment variables in the LiveOcean system
    using values in the dict "Ldir".  It is similar to Lstart.m in the
    MATLAB code, but it returns a dictionary instead of a structure.

    We use input parameters to allow for different gridnames and tags to
    coexist.
    """
    # put top level information from input into a dict
    Ldir = dict()
    Ldir['gridname'] = gridname
    Ldir['tag'] = tag
    if os.path.isfile(alp_user + '/get_lo_info.sh'):
        subprocess.call([alp_user + '/get_lo_info.sh'])
    else:
        subprocess.call([alp + '/get_lo_info.sh'])
    Ldir_temp = csv_to_dict(alp + '/lo_info.csv')
    Ldir.update(Ldir_temp)
    # and add a few more things
    Ldir['gtag'] = Ldir['gridname'] + '_' + Ldir['tag']
    Ldir['grid'] = Ldir['data'] + 'grids/' + Ldir['gridname'] + '/'
    Ldir['forecast_days'] = 3
    Ldir['LOu'] = upath + '/'
    return Ldir

def make_dir(dirname, clean=False):
    # Make a directory if it does not exist.
    # Use clean=True to clobber the existing directory.
    if clean == True:
        import shutil
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass # assume OSError was raised because directory already exists

def dict_to_csv(dict_name, csv_name, write_mode='w'):
    # Write the contents of a dict to a two-column csv file.
    # The write_mode can be w (overwrite) or a (append).
    with open(csv_name, write_mode) as ff:
        ww = csv.writer(ff)
        for kk in dict_name.keys():
            ww.writerow((kk, dict_name[kk]))

def csv_to_dict(csv_name):
    # Reads two-column csv file into a dict.
    dict_name = dict()
    with open(csv_name) as ff:
        for row in csv.reader(ff):
            dict_name[row[0]] = row[1]
    return dict_name

def run_worker(Ldir, worker_type='matlab'):
    # run the worker code using subprocess
    if worker_type == 'matlab':
        # pass arguments to a matlab program
        func = ("make_forcing_worker(\'" +
            Ldir['gridname'] + "\',\'" +
            Ldir['tag'] + "\',\'" +
            Ldir['date_string'] + "\',\'" +
            Ldir['run_type'] + "\',\'" +
            Ldir['LOogf_f'] + "\')")
        cmd = Ldir['which_matlab']
        run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
        proc = subprocess.run(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('\n-main: screen output from worker-')
        print(proc.stdout.decode())
    else:
        print('other worker types not implemented yet')

def datetime_to_modtime(dt):
    # This is where we define how time will be treated
    # in all the model forcing files.
    # NOTE: still need to define time zone (UTC?)

    # INPUT:
    # dt is a single datetime value

    # this returns seconds since 1/1/1970
    t = (dt - datetime(1970,1,1,0,0)).total_seconds()
    return t

def modtime_to_datetime(t):
    # INPUT: seconds since 1/1/1970 (single number)
    # OUTPUT: datetime version
    dt = datetime(2008,1,1,0,0) + timedelta(seconds=t)
    return dt

def modtime_to_mdate_vec(mt_vec):
    # INPUT: numpy vector of seconds since 1/1/1970
    # mt stands for model time
    # OUTPUT: a vector of mdates
    # first make a list of datetimes
    import matplotlib.dates as mdates
    dt_list = []
    for mt in mt_vec:
        dt_list.append(datetime(1970,1,1,0,0) + timedelta(seconds=mt))
    md_vec = mdates.date2num(dt_list)
    return md_vec

# Functions used by postprocessing code like pan_plot or the various extractors

def date_list_utility(dt0, dt1):
    # INPUT: start and end datetimes
    # OUTPUT: list of LiveOcean formatted dates
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + timedelta(1)
    return date_list

def fn_list_utility(dt0, dt1, Ldir, hourmax=24):
    # INPUT: start and end datetimes
    # OUTPUT: list of all history files expected to span the dates
    dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
    fn_list = []
    date_list = date_list_utility(dt0, dt1)
    for dl in date_list:
        f_string = 'f' + dl
        if dl == date_list[0]:
            hourmin = 0
        else:
            hourmin = 1
            # skip hour zero on subsequent days
            # because it is a repeat
        for nhis in range(hourmin+1, hourmax+2):
            nhiss = ('0000' + str(nhis))[-4:]
            fn = dir0 + f_string + '/ocean_his_' + nhiss + '.nc'
            fn_list.append(fn)
    return fn_list

def get_fn_list(list_type, Ldir, date_string0, date_string1, his_num=1):
    dt0 = datetime.strptime(date_string0, '%Y.%m.%d')
    dt1 = datetime.strptime(date_string1, '%Y.%m.%d')
    dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
    if list_type == 'snapshot':
        # a single file name in a list
        his_string = ('0000' + str(his_num))[-4:]
        fn_list = [dir0 + 'f' + date_string0 +
                   '/ocean_his_' + his_string + '.nc']
    elif list_type == 'hourly':
        # list of hourly files over a date range
        fn_list = fn_list_utility(dt0,dt1,Ldir)
    elif list_type == 'forecast':
        # list of all history files in a directory
        hourmax = Ldir['forecast_days'] * 24
        dt0 = datetime.now()
        fn_list = fn_list_utility(dt0, dt0, Ldir, hourmax=hourmax)
    elif list_type == 'low_passed':
        # list of low passed files over a date range
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = (dir0 + f_string + '/low_passed.nc')
            fn_list.append(fn)
    elif list_type == 'daily':
        # list of history file 21 over a date range
        # which is Noon Pacific Standard Time
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = (dir0 + f_string + '/ocean_his_0021.nc')
            fn_list.append(fn)
    elif list_type == 'daily4':
        # list of 4 history files per day over a date range
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = (dir0 + f_string + '/ocean_his_0003.nc')
            fn_list.append(fn)
            fn = (dir0 + f_string + '/ocean_his_0009.nc')
            fn_list.append(fn)
            fn = (dir0 + f_string + '/ocean_his_0015.nc')
            fn_list.append(fn)
            fn = (dir0 + f_string + '/ocean_his_0021.nc')
            fn_list.append(fn)
    elif list_type == 'merhab':
        # a list of all but the first history file in a directory
        in_dir = dir0 + 'f' + date_string0 + '/'
        fn_list_raw = os.listdir(in_dir)
        fn_list = [(in_dir + ff) for ff in fn_list_raw if 'ocean_his' in ff]
        fn_list.sort()
        fn_list.pop(0) # remove the first hour
    elif list_type == 'allhours':
        # a list of all the history files in a directory
        in_dir = dir0 + 'f' + date_string0 + '/'
        fn_list_raw = os.listdir(in_dir)
        fn_list = [(in_dir + ff) for ff in fn_list_raw if 'ocean_his' in ff]
        fn_list.sort()

    return fn_list

def choose_item(indir, tag='', exclude_tag='', itext='** Choose item from list **'):
    print('\n%s\n' % (itext))
    ilist_raw = os.listdir(indir)
    ilist_raw.sort()
    if len(tag) == 0:
        ilist = [item for item in ilist_raw if item[0] != '.']
    else:
        ilist = [item for item in ilist_raw if tag in item]

    if len(exclude_tag) == 0:
        pass
    else:
        ilist = [item for item in ilist if exclude_tag not in item]

    Nitem = len(ilist)
    idict = dict(zip(range(Nitem), ilist))
    for ii in range(Nitem):
        print(str(ii) + ': ' + ilist[ii])
    my_choice = input('-- Input number -- (return=0)')
    if len(my_choice)==0:
        my_choice = 0
    my_item = idict[int(my_choice)]
    return my_item

def copy_to_azure(input_filename, output_filename, container_name, Ldir):
    """
    Copy the file to azure, and return the URL to access the file.
    Use the new azure module structure from "pip install azure-storage-blob".
    2020_06
    """
    # get the blob_service
    from azure.storage.blob import BlobServiceClient
    from azure.storage.blob import PublicAccess
    azu_dict = csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
    account = azu_dict['account']
    key = azu_dict['key']
    connection_string = ('DefaultEndpointsProtocol=https' +
        ';AccountName=' + account +
        ';AccountKey=' + key +
        ';EndpointSuffix=core.windows.net')
    blob_service = BlobServiceClient.from_connection_string(conn_str=connection_string)
    # create the container if needed
    try:
        blob_service.create_container(container_name, public_access=PublicAccess.Container)
    except:
        # assume error is because container exists
        pass
    # write it to Azure and return a dict of information
    az_dict = {}
    try:
        from azure.storage.blob import BlobClient
        blob = BlobClient.from_connection_string(conn_str=connection_string,
            container_name=container_name, blob_name=output_filename)
        with open(input_filename, 'rb') as data:
            blob.upload_blob(data, overwrite=True)
        az_dict['result'] = 'success'
        az_url = ('https://pm2.blob.core.windows.net/'
            + container_name + '/' + output_filename)
        az_dict['az_url'] = az_url # URL to get the file
    except Exception as e:
        az_dict['result'] = 'fail'
        az_dict['exception'] = str(e)
    return az_dict
