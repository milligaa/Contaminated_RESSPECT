import numpy as np
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import glob
from scipy import integrate, interpolate, ndimage
from qmostetc import Spectrum, QMostObservatory, L1DXU



def setup_data(data):

    #extract all the data we need, this will now include the ddlr and snsep values required to calculate the effective fibre magnitude

    Smags = []
    Gmags = []
    snid = []
    ddlr = []
    snsep = []
    texp_visit = []

    real_data = data

    for it in range(len(real_data)):
        Smags.append(real_data['sne_rband'][it])
        Gmags.append(real_data['gal_rband'][it])
        snid.append(real_data['snid'][it])
        ddlr.append(real_data['ddlr'][it])
        snsep.append(real_data['snsep'][it])
        texp_visit.append(real_data['texp_seconds'][it] / 60)



    galaxies = []
    supernovae = []



    for t in snid:
        file_location = '/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/RESSPECT/Spectra_120224/cid'+str(t)
        both_spectrum = glob.glob(str(file_location)+'/*.SPEC')
        host_spectrum = glob.glob(str(file_location)+'/*HOST.SPEC')
        trans_spectrum = []

        for file in both_spectrum:
            if "HOST" in file:
                continue
            else:
                trans_spectrum.append(file)

        supernovae.append(trans_spectrum[0])
        galaxies.append(host_spectrum[0])

    return Smags, Gmags, galaxies, supernovae, ddlr, snsep, texp_visit, snid

def host_fibremag(sep, gal_ddlr, sersic_index, galmag, pix_size_arcsec, seeing):

    dlr = gal_ddlr * sep
    gmag = galmag
    seperation = sep
    
    #first get total galaxy flux in erg/s/cm2 (constant here is the zero-point flux in AB system in the LSST r-band filter bandpass)
    flux_gal_total = 4.69542e-6 * (np.e ** (-gmag / 2.5))

    #this flux is equal to 2.8941 * pi * Ie * Re^2 (but only in the case of sersic index of 0.5)
    eff_intensity = flux_gal_total / (2.8941 * np.pi * (dlr ** 2))

    #then calculate effective intensity at centre of fibre (seperation)
    bn = (2 * sersic_index) - (1/3) + (0.009876 * sersic_index)  #from Prugneil and Simien

    #now create the pixel array from integer number of pixels in display range
    pixel_length = pix_size_arcsec
    pixel_no = int((6 * seperation) / pixel_length) + 1 #plus one because always rounds down

    int_array = []
    distance_array = []
    for i in range(pixel_no):
        sub_array = []
        sub_dist = []
        for j in range(pixel_no):
            sub_array.append(0)
            sub_dist.append(np.sqrt((i - (pixel_no/2))**2 + (j - (pixel_no/2))**2) * pixel_length)
        
        int_array.append(sub_array)
        distance_array.append(sub_dist)

    #define the intensities at each pixel
    for h in range(pixel_no):
        for g in range(pixel_no):
            int_array[h][g] = np.log(eff_intensity * np.e ** (-bn * (((distance_array[h][g] / dlr) ** (1/sersic_index)) - 1)))

    #define min and max intensities used in iamge (not max and min overall) and other useful numbers, objects
    min_int = min(min(int_array))

    #define conv radius using seeing (FWHM -> sigma)
    conv_radius = (seeing / np.sqrt(8 * np.log(2))) / pixel_length
    new_int_array = ndimage.gaussian_filter(int_array, sigma=conv_radius, mode='constant', cval=min_int)

    pixel_sep = seperation / pixel_length
    pixel_fibre_size = 0.725 / pixel_length
    true_array_int = np.asarray(new_int_array)

    #now try and get the pixels that are in the fibre. pixel coords are to the bottom left corner, but centres are at integer coords
    #remember circle of fibre is N/2 pixels up and to right 
    good_coords = []
    fibre_int_pix = 0
    for h in range(pixel_no):
        for g in range(pixel_no):
            pixel_centre_x = h
            pixel_centre_y = g
            pixel_distance = (pixel_centre_x - ((pixel_no/2) + pixel_sep)) ** 2 + (pixel_centre_y - (pixel_no/2)) ** 2
            pixel_int = true_array_int[h][g]
            if pixel_distance <= pixel_fibre_size ** 2:
                # print((np.e ** pixel_int) * pixel_length * pixel_length, 'this is the pixel flux for pixel', h, g)
                fibre_int_pix = fibre_int_pix + ((np.e ** pixel_int) * pixel_length * pixel_length)
                good_coords.append([h,g])
            else:
                continue

    for i in range(len(good_coords)):
        coords = good_coords[i]
        int_array[coords[1]][coords[0]] = min_int


    ratio_pix = fibre_int_pix / flux_gal_total

    if ratio_pix > 1:
        ratio_pix = 1
    else:
        ratio_pix = ratio_pix

    #now calculate effective galaxy mag
    eff_mag = gmag - 2.5 * np.log10(ratio_pix)

    return eff_mag


def transient_fibremag(seeing, sne_mag):
    #first turn the seeing value into a sigma for the gaussian
    FWHM = seeing
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))

    #generate a radial guassian (gaussian multiplied by 2 pi x)
    gaussian = lambda x:2 * np.pi * x * np.e ** (-(x**2)/(2 * (sigma ** 2)))

    #integrate to infinity to get the normalisation
    normalisation = integrate.quad(gaussian, 0, np.inf)

    #integrate again to fibre radius with normalisation constant dividing gaussian
    gaussian_norm = lambda x:2 * np.pi * x * np.e ** (-(x**2)/(2 * (sigma ** 2))) / normalisation[0]
    fraction = integrate.quad(gaussian_norm ,0 ,0.725)[0] #only want the first bit here

    print(fraction)

    new_mag = sne_mag - 2.5*np.log(fraction)
    return new_mag


def ETC_specMaker(SNe_data, Gal_data, Gal_mag, SNe_mag, SN_type, texp, seeing, spec_save_path):
    #this program calculates the mag of data, then sets a template equal to that mag
    #then it interpolates the template flux over the data's wavelength values before contaminating the data with the template flux
    #plots can be made of the spectra at each stage in the pipeline and TiDES SNR of the raw output can be calculated
    #final L1 combined spectrum is saved (wavelength, flux, sky flux, error flux and sky + error flux I think (don't quote me on that (i don't know what I'm doing)))

    #defining parameters for use down the code
    G_mag = Gal_mag * u.ABmag
    S_mag = SNe_mag * u.ABmag


    #importing the data to be used

    f = open(SNe_data, "r")
    f_lines = f.readlines()
    lam = []
    intens = []
    for i in range(len(f_lines)):
        if not f_lines[i].startswith("#"):
            ls = f_lines[i].split()
            lam.append(float(ls[0]))
            intens.append(float(ls[1]))

    #making a pectrum object and then using the magnitude of the sectrum to set a new magnitude

    wavelength = lam * u.nm / 10
    photon_flux = intens * u.erg / (u.cm ** 2 * u.s * u.angstrom)
    wavelength_z = wavelength

    spec = Spectrum(wavelength_z, photon_flux)
    spec_mag = spec.get_mag(u.ABmag, 'LSST_LSST.r')

    #note the 0.545 mag factor added to Smag to account for fibre loss
    flux_ratio = 10 ** ((S_mag.value - spec_mag.value)/(-2.5)) 

    new_flux = photon_flux * flux_ratio

    #okie dokie now lets try and repeat the above for the template data

    f2 = open(Gal_data, "r")
    f2_lines = f2.readlines()
    lam2 = []
    intens2 = []
    for i in range(len(f2_lines)):
        if not f2_lines[i].startswith("#"):
            ls2 = f2_lines[i].split()
            lam2.append(float(ls2[0]))
            intens2.append(float(ls2[1]))

    #only take a part of the spectrum, too wide a range breaks the interpolation later
    wavelength2 = lam2 * u.nm / 10     #divide by ten to 'make' wavelength nm and agree with units of SNe data
    photon_flux2 = intens2 * u.erg / (u.cm ** 2 * u.s * u.angstrom)
    wavelength_z2 = wavelength2

    spec2 = Spectrum(wavelength_z2, photon_flux2)
    spec_mag2 = spec2.get_mag(u.ABmag, 'LSST_LSST.r')

    flux_ratio2 = 10 ** ((G_mag.value - spec_mag2.value)/(-2.5)) 

    new_flux2 = photon_flux2 * flux_ratio2

    #now we interpolate the galaxy spectrum on to the same wavelength values as the SNe spectrum
    #we start by defining the redshifted spectrum eavelengths here

    g_wave_z = spec2.wavelength
    g_flux_value = new_flux2.value
    g_wave_value = g_wave_z.value        

    #now we perform the interpolation

    interp = interpolate.splrep(g_wave_value, g_flux_value)
    interp_gal = interpolate.splev(wavelength_z.value, interp)


    #now we add varing amounts of contamination from galaxy to SN
    contaminated_flux = []

    for i in range(len(lam)):
        interp_gal[i] = interp_gal[i]
        contaminated_flux.append(new_flux[i].value + interp_gal[i])


    
    # Object to simulate the 4MOST observatory, including atmosphere,
    # telescope, spectrograph, CCD.
    qmost = QMostObservatory('lrs')
    obs = qmost(33.55730976*u.deg, seeing*u.arcsec, 'dark')

    #define a new spectrum obk=ject with the redshifted wavelength and the SN flux with host galaxy contamination and set up correct spectral units
    new_spec = Spectrum(wavelength_z, contaminated_flux)
    new_spec.flux= new_spec.flux * u.erg / (u.cm ** 2 * u.s * u.angstrom)
    new_mag = new_spec.get_mag(u.ABmag, 'LSST_LSST.r')
    #divide by arcseconds squared for flat spatial distro, already account for fibre size by adding 0.545 mag to SNe
    new_spec.flux = new_spec.flux / (1.665 * u.arcsec * u.arcsec)

    

    #set the new spectrum as the target for observation and then set the observational parameter

    obs.set_target(new_spec, 'flat')
    res = obs.expose(texp*u.min, 1) 

    #hopefully this little bit of code here will produce the out put but with all 3 arms joined
    #we create a wrapper for the output (called dxu here) and then write the joined data to a fits file

    dxu = L1DXU(qmost, res, texp*u.min)   #doesn't like exposures sperarated into blocks, not sure if matters
    comb_spec = dxu.joined_spectrum()

        
    #first import the transmission spectrum
    T_data = '/Users/andrew/Desktop/Python_Stuff/SN_and_Galaxy/adjusted_flat_spec.txt'
    T_data_table = Table.read(T_data, format = 'csv', delimiter = ' ')

    lam_t = []
    intens_t = []
    for rty in range(len(T_data_table)):
        lam_t.append(T_data_table[rty][0])
        intens_t.append(T_data_table[rty][1])

    #for now the T spectrum is generated using a combined spectrum so wavelength array matches without need for interpolation
    #THIS MAY NOT ALWAYS BE TRUE, BE WARY TRAVELLER
    corrected_flux = []
    for q in range(len(lam_t)):
        corrected_flux.append(comb_spec['FLUX'][q].value / intens_t[q])

    comb_spec['FLUX'][:] = corrected_flux * u.erg / (u.cm ** 2 * u.s * u.angstrom)

    #setting up strings for data wanted in the file names for writing the combined spec
    st = str(SN_type)
    sm = str(round(SNe_mag, 5))
    gm = str(round(Gal_mag, 5))
    exp = str(round(texp, 5))
    
    #exreacting the data from the combined spectra

    comb_table = Table()
    comb_table['wavelength'] = comb_spec['WAVE']
    comb_table['flux'] = comb_spec['FLUX']
    comb_table['err'] = comb_spec['ERR_FLUX']
    comb_table['flux_nos'] = comb_spec['FLUX_NOSS']
    comb_table['err_nos'] = comb_spec['ERR_FLUX_NOSS']

    #setting up the system path for saving the file in the right place and then writing it there
    
    name_of_file = spec_save_path+st+'Smag'+sm+'Gmag'+gm+'texp'+exp+'.txt'
    
    ascii.write(comb_table, name_of_file, format = 'csv', overwrite = False)
    
    #calculating the combined galaxy and SNe SNR for the contaminated spectrum
    comb_SNR_unbin = []
    for no in range(len(comb_spec['WAVE'])):
        comb_SNR_unbin.append(comb_spec['FLUX'][no]/comb_spec['ERR_FLUX_NOSS'][no]) 


    #print(comb_spec['WAVE'])
    snr_bin = []
    noise_bin = []
    target_bin = []
    wl_bin = []
    bin_number = int((17315 - 3315)/ 60)
    #print('this is bin number', bin_number)
    for ijk in range(bin_number):
        noise_sum = 0
        target_sum = 0
        wl_sum = 0
        tot = 0
        for jki in range(60):
            target_sum = target_sum + comb_spec['FLUX'][3315 + (ijk * 60) + jki].value
            noise_sum = noise_sum + (comb_spec['ERR_FLUX_NOSS'][3315 + (ijk * 60) + jki].value ** 2)
            wl_sum = wl_sum + comb_spec['WAVE'][3315 + (ijk * 60) + jki].value
            tot = tot + 1

        snr_bin.append(target_sum / np.sqrt(noise_sum))
        noise_bin.append(np.sqrt(noise_sum) / tot)
        target_bin.append(target_sum / tot)
        wl_bin.append(wl_sum / tot)

    comb_snr = np.sum(snr_bin) / bin_number

    return [comb_snr, new_mag, comb_spec['WAVE'], comb_spec['FLUX']]