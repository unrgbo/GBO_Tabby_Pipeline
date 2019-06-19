import Tabby_pipe as tpipe
from tqdm import tqdm
import datetime

print '\nThe current datetime is:\n'
print str(datetime.datetime.now())

print '\nLooking for new data from the GBO\n'

tpipe.move_new_files()

donedates, baddates, obsdates, biasdates, darkdates, flatdates = tpipe.get_dates(band='r')

newdates = [x for x in obsdates if x not in donedates]
newdates = [x for x in newdates if x not in baddates]


if len(newdates) == 0:
    print '\nDid not find any new data to process\n\nExiting Tabby_API'

elif len(newdates) >0:
    print '\nBeginning data processing for new dates\n'

    for date in tqdm(newdates):

        print '\nBeginning photometry for {}'.format(date)
        tpipe.do_phot(date)

        print '\nParsing output from IRAF photometry files\n'
        tpipe.night_phot(date)

        print '\nAdding data from {} to cumulative data file\n'.format(date)
        tpipe.add_new_dates()

        print '\nMaking summary of all Tabby data and output\n'
        tpipe.make_summary()

        print '\nMaking updated plot and writing to pipe_out\n'
        tpipe.make_plot(write = True)

    print '\nBeginning backups\n'
    tpipe.backup()
