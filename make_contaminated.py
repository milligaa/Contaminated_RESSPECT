#just a version of the comb_spec maker that produces SNR values for missing objects
from astropy.io import ascii
from astropy.table import Table
import sys
import os
from datetime import date

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from util.extract_RESSPECT_data import get_RESSPECT_data
from util.utils import setup_data, host_fibremag, transient_fibremag, ETC_specMaker

def template_flow(input_path, output_path):

    #generates params.csv file from RESSPECT sepctra

    get_RESSPECT_data(input_path)

    #uses the params.csv to get parameters like redshift, magnitude and SNe type for use in template generation
    template_values = setup_data(Table.read(input_path+'params.csv', format='csv', delimiter=','), input_path)
    Smags = template_values[0]
    Gmags = template_values[1]
    galaxies = template_values[2]
    supernovae = template_values[3]
    ddlr = template_values[4]
    snsep = template_values[5]
    texp_visit = template_values[6]
    snid = template_values[7]

    #now must perform the correction for effective fibre mag
    seeing_val = 0.8
    gmag_eff_fibre = []
    smag_eff_fibre = []
    for h in range(len(Gmags)):
        gmag_eff_fibre.append(host_fibremag(snsep[h], ddlr[h], 0.5, Gmags[h], snsep[h]/100, seeing_val))
        smag_eff_fibre.append(transient_fibremag(seeing=seeing_val, sne_mag=Smags[h]))

    #set the host and transient magnitudes to be the calculated fibre magnitudes
    for o in range(len(gmag_eff_fibre)):
        Gmags[o] = gmag_eff_fibre[o]
        Smags[o] = smag_eff_fibre[o]


    L1_SNR_corr1 = []
    Comb_mag_corr1 = []

    bad_index = []
    SNR_append_index = []

    SN_type_str = []

    for j in range(len(Gmags)):
        SN_type_str.append(snid[j])

    for dummy in range(len(Gmags)): 
            
            print(dummy + 1, 'of', len(supernovae))
            
            result_observed = ETC_specMaker(supernovae[dummy], galaxies[dummy],
            Gmags[dummy], Smags[dummy], SN_type_str[dummy], texp_visit[dummy],
            seeing_val, output_path)
    
            L1_SNR_corr1.append(result_observed[0])
            Comb_mag_corr1.append(result_observed[1].value)

            SNR_append_index.append(dummy)


    for bad in range(len(bad_index)):
        del(Smags[bad_index[bad] - bad])
        del(Gmags[bad_index[bad] - bad])
        del(texp_visit[bad_index[bad] - bad])

    SNR_table = Table()
    SNR_table['SN_ID'] = snid
    SNR_table['combined_SNR'] = L1_SNR_corr1
    SNR_table['combined_mag'] = Comb_mag_corr1
    SNR_table['fibre_trans_mag'] = Smags
    SNR_table['fibre_host_mag'] = Gmags
    SNR_table['esposure'] = texp_visit

    #saves some output parameters for the generated spectra in a csv file with
    #the current date (ddmmYYYY) in the filename
    today = date.today()

    d1 = today.strftime("%d%m%Y")
    ascii.write(SNR_table, output_path+'output_params'+d1+'.csv', format='csv', overwrite = True)

if __name__ == "__main__":
    template_flow('input_bank/', 'output_bank/')