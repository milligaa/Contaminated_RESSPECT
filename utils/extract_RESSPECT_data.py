# looks into folder in the format of Spectra_120224 and gets all the information
# required to pass onto RESSPECT_maker.py

import glob
from astropy.table import Table
from astropy.io import ascii

def get_RESSPECT_data(path):
    SNe_events = glob.glob(str(path)+'/*')

    for i in range(len(SNe_events)):
        if '.csv' in SNe_events[i]:
            del(SNe_events[i])
        else:
            continue

    snsep = []
    ddlr = []
    gal_mag = []
    sne_mag = []
    snid = []
    texp = []

    for folder in SNe_events:
        both_spectrum = glob.glob(str(folder)+'/*.SPEC')
        host_spectrum = glob.glob(str(folder)+'/*HOST.SPEC')
        data_file = glob.glob(str(folder)+'/*.DAT')
        trans_spectrum = []

        print(folder)

        for file in both_spectrum:
            if "HOST" in file:
                continue
            else:
                trans_spectrum.append(file)
        
        snid_b = data_file[0].find('SALT2_cid') + 9
        snid_e = data_file[0].find('_mjd')
        snid.append(data_file[0][snid_b:snid_e])

        f = open(data_file[0], "r")
        data_text = f.readlines()
        for i in range(len(data_text)):
            if 'HOSTGAL_SNSEP' in data_text[i]:
                ls = data_text[i].split()
                snsep.append(ls[1])

            elif 'HOSTGAL_DDLR' in data_text[i]:
                ls = data_text[i].split()
                ddlr.append(ls[1])

            elif 'HOSTGAL_MAG:' in data_text[i]:
                ls = data_text[i].split()
                gal_mag.append(ls[3])

            elif 'SPECTRUM_MJD' in data_text[i]:
                if data_text[i].split()[3] == 'HOST':
                    continue
                else: 
                    sne_mag.append(data_text[i+6].split()[3])

            elif 'SPECTRUM_ID' in data_text[i]:
                if int(data_text[i].split()[1]) == 2:
                    texp.append(data_text[i + 2].split()[1])
                else:
                    continue

    print(len(texp), len(snid))

    to_save = Table()
    to_save['snid'] = snid
    to_save['ddlr'] = ddlr
    to_save['snsep'] = snsep
    to_save['gal_rband'] = gal_mag
    to_save['sne_rband'] = sne_mag
    to_save['texp_seconds'] = texp


    ascii.write(to_save, path+'/params.csv', format = 'csv', overwrite=True)