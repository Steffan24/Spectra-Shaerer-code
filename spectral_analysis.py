

from variables import output_file, output_flux, output_std, output_scale, dir_basic, output_mass, output_exp_time,  mass_output_array, input_scale, dir_basic_logA, dir_basic_logE, dir_basic_high_res
from modules import np, plt, curve_fit, cm, ascii, cosmo
import plotting_params
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

output_exp_time = 360
output_mass = 4.5

spectrum = np.load(f"{dir_basic_logE}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_counts.npy")
spectrum_flux = np.load(f"{dir_basic_logE}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_flux.npy")
spectrum_std = np.load(f"{dir_basic_logE}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_std.npy")
wavelength_angstrom = np.load(f"{dir_basic_logE}/{output_exp_time}ks_exposures/M_{output_mass}_output/wavelength_angstrom.npy")

def ks_to_seconds(t_ks):
    return t_ks * 1000.0

def gaussian(x,A):
    std = 3.244
    B = 20800
    gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
    return gaussian

def gaussian_line(x,A, B, C):
    std = C
    B = B
    gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
    return gaussian

def chi_squared_gauss(x,y, std, popt):
    chi = ((y - gaussian_line(x, *popt))**2)/(std**2)
    return np.sum(chi)

def chi_squared_cont(x,y,std, straight_fit):
    chi = ((y - straight_fit)**2) / (std**2)
    #print(chi)
    return np.nansum(chi)

def calcSNR_chi2(wavelength_angstrom,spectrum, spectrum_std):
    mask_gauss = (wavelength_angstrom > 19000) & (wavelength_angstrom < 22000)
    wavelength_gauss = wavelength_angstrom[mask_gauss]
    counts_gauss = spectrum[mask_gauss]
    std_gauss = np.array(spectrum_std[mask_gauss])
    # plt.figure()
    # plt.step(wavelength_gauss, counts_gauss)
    # plt.show()
    print(spectrum_std.shape)
    center = 20800
    sigma_guess = 3.244
    p0 = [np.nanmax(counts_gauss) - np.nanmedian(counts_gauss), center, sigma_guess]  # [amp, mu, sigma]
    bounds = ([0, center-1.5, 3], [np.inf, center+1.5, 4])
    popt, pcov, infodict, mesg, ier = curve_fit(gaussian_line, wavelength_gauss,counts_gauss,sigma = std_gauss,absolute_sigma = True, p0=p0, bounds=bounds, full_output = True)
    #print(f"std gauss: {popt[1]}")
    print(f"popt: {popt}")
    print(mesg)
    print(ier)
    chi_squared_gaussian = chi_squared_gauss(wavelength_gauss,counts_gauss, std_gauss, popt)
    print(f"chi squared gauss: {chi_squared_gaussian}")

    mask_line = (wavelength_gauss > 20804) | (wavelength_gauss < 20796)
    wavelength_line = wavelength_gauss[mask_line]
    counts_line = counts_gauss[mask_line]
    std_line = std_gauss[mask_line]
    straight_fit = np.ones(len(wavelength_gauss)) * np.nanmedian(np.array(counts_line))
    #print(f"straight fit: {straight_fit}")
    chi_squared_cont_value = chi_squared_cont(wavelength_gauss, counts_gauss, std_gauss, straight_fit)
    print(f"chi squared cont: {chi_squared_cont_value}")

    SNR = np.sqrt(np.abs(chi_squared_gaussian - chi_squared_cont_value))
    print(f"SNR_chi_method: {SNR}")
    total_counts = popt[0]
    
    # plt.plot(wavelength_gauss, counts_gauss)
    # plt.show()
    # plt.figure()
    # plt.plot(wavelength_gauss, counts_gauss)
    # plt.plot(wavelength_gauss, straight_fit)
    # plt.show()
    # plt.figure()
    # plt.plot(wavelength_angstrom, spectrum)
    # plt.show()
    narrow = (wavelength_gauss < 20850) & (wavelength_gauss > 20750)
    narrow_wave = wavelength_gauss[narrow]
    counts_narrow = counts_gauss[narrow]
    # plt.figure(111)
    # plt.plot(narrow_wave, gaussian_line(narrow_wave, *popt), c = 'red', label='$Gaussian\ fit$')
    # plt.step(narrow_wave, counts_narrow)
    # plt.xlim(20750, 20850)
    # plt.xlabel(rf"$Wavelength\ (\mathring{{A}}) $")
    # plt.ylabel(rf'$Counts$')
    # plt.legend(loc='upper left', bbox_to_anchor=[0.1, 0.9, 0.4, 0.3])
    # plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/gaussian.png')
    # plt.show()
    return SNR, total_counts


def calcSNR_peak(wavelength_angstrom, spectrum, spectrum_std):
    mask_gauss =  (wavelength_angstrom > 20797) & (wavelength_angstrom < 20803)
    pixel_sizes = np.diff(wavelength_angstrom)
    print("Median pixel scale:", np.median(pixel_sizes))
    print("Unique values:", np.unique(np.round(pixel_sizes, 3)))
    wavelength_gauss = wavelength_angstrom[mask_gauss]
    counts_gauss = spectrum[mask_gauss]
    std_gauss = spectrum_std[mask_gauss]
    # plt.figure()
    # plt.step(wavelength_gauss, counts_gauss)
    # plt.show()
    # popt, pcov, infodict, mesg, ier = curve_fit(gaussian, wavelength_gauss,counts_gauss,sigma = std_gauss,absolute_sigma = True, p0=[1*10**5], full_output = True)
    # A = popt[0]
    # peak = A / (np.sqrt(2*np.pi)*(3.244))
    # mask_cont = (wavelength_angstrom > 18000) & (wavelength_angstrom < 23000)
    # spectrum_cont = spectrum[mask_cont]
    # wavelength_cont = wavelength_angstrom[mask_cont]
    # std = np.nanstd(spectrum_cont[(wavelength_cont < 19500) | (wavelength_cont > 21000)])
    # SNR = peak / std*np.sqrt(2.355*(3.244/2.642))
    # total_counts = A

    sum_counts = np.nansum(counts_gauss)
    sum_std = np.sqrt(np.nansum(std_gauss**2))
    SNR = sum_counts/sum_std
    print(f"chi squared method 2: {SNR}")
    return SNR

print(f"SNR method 2: {calcSNR_peak(wavelength_angstrom, spectrum, spectrum_std)}")

def calc_line_flux(spectrum_flux, spectrum_std, wavelength_angstrom):
    mask_gauss =  (wavelength_angstrom > 20795) & (wavelength_angstrom < 20805)
    wavelength_gauss = wavelength_angstrom[mask_gauss]
    flux_gauss = spectrum_flux[mask_gauss]
    std_gauss = np.array(spectrum_std[mask_gauss])
    print(spectrum_std.shape)
    popt, pcov, infodict, mesg, ier = curve_fit(gaussian_line, wavelength_gauss,flux_gauss,sigma = std_gauss, p0=[1*10**(-17), 20800, 3.244], full_output = True)
    line_flux = popt[0]
    flux_err = np.sqrt(np.abs(pcov[0, 0]))
    print("popt =", popt)
    print("pcov =", pcov)
    print("sqrt diag =", np.sqrt(np.diag(pcov)))
    # plt.figure()
    # plt.plot(wavelength_gauss, flux_gauss)
    # plt.plot(wavelength_gauss, gaussian(wavelength_gauss, *popt))
    # plt.show()
    
    return line_flux, flux_err

#def calc_theory_line_flux(mass_output_array):
    
    


def collapse_all(mass_output_array, output_exp_time):
    n_files = len(mass_output_array)
    fig, axes = plt.subplots(n_files, 2, figsize=(10, 3 * n_files),
                             gridspec_kw={'width_ratios': [1.5, 1]})
    #axes = axes.flatten()
    for i in range(len(mass_output_array)):
        data = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/data.npy")
        spectrum = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/spectrum_counts.npy")
        flux = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/spectrum_flux.npy")
        spectrum_std = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/spectrum_std.npy")
        data_std = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/data_std.npy")
        data_flux = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/data_flux.npy")
        wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass_output_array[i]}_output/wavelength_angstrom.npy")

        # plt.figure()
        # plt.plot(wavelength_angstrom, spectrum)
        # plt.show()
        
        zoom_center = 20800
        zoom_halfwidth = 50
        mask = (wavelength_angstrom > zoom_center - zoom_halfwidth) & (wavelength_angstrom < zoom_center + zoom_halfwidth)
        collapsed_image = np.nansum(data[mask,:,:],axis=0)
        collapsed_image = np.nan_to_num(collapsed_image)
        spectrum_cut = spectrum[mask]
        wavelength_cut = wavelength_angstrom[mask]

        zoom  = axes[i,0]
        zoom.plot(wavelength_cut, spectrum_cut*10**(-6))
        popt, pcov = curve_fit(gaussian, wavelength_cut,spectrum_cut, p0=[1])
        zoom.plot(wavelength_cut, gaussian(wavelength_cut, *popt)*10**(-6), ls='--', c='red')
        zoom.set_xlim(zoom_center - zoom_halfwidth, zoom_center + zoom_halfwidth)
        if i == n_files - 1:
            zoom.set_xlabel("\(\lambda (\mathring{A})\)")
        else:
            zoom.set_xticklabels([])
            #zoom.set_yticks([])
            #zoom.set_xticks([])
        
        #zoom.set_ylim(-0.01, 0.01)
        im2_ax = axes[i,1]
        ny, nx = collapsed_image.shape
        pixel_scale_mas = 7
        extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
        im2 = im2_ax.imshow(collapsed_image, origin='lower', cmap='gray',
                            aspect='equal', extent=extent)
        constants = [8.0,7.5, 7, 6.5,6.0, 5.5,5.0,4.5,4.0]
        constants = np.array(constants)
        if i == n_files - 1:
            im2_ax.set_xlabel("\(X (mas)\)")
        else:
            im2_ax.set_xticklabels([])
            im2_ax.set_yticklabels([])
            #im2_ax.set_xticks([])
        im2_ax.axis('on')

        flux_cut = flux[mask]
        F_array = np.median(flux_cut)
        wavelength_ref = np.median(wavelength_cut)
        V_median = -2.5*np.log10(F_array) - 2.402 - 5.0*np.log10(wavelength_ref)
        SNR_val, total_counts = calcSNR_chi2(wavelength_angstrom,spectrum, spectrum_std)
        im2_ax.text(1.2, 0.5,fr"$M_{{cluster}} = 10^{{{mass_output_array[i]:.2f}}} M_{{\odot}}$"
                    "\n"
                    fr"$SNR = {SNR_val:.2f}$",
                    transform=im2_ax.transAxes, ha='left', va='center', fontsize=16)
       # im2_ax.text(1.2, 0.5,fr"$M_{{cluster}} = 10^{{{constants[i]:.2f}}} M_{{\odot}}$"
#"\n"
#fr"$V_{{median}} = {V_median:.2f}$"
 #                   "\n"
  #                  fr"$SNR = {SNR_val:.2f}$",
   #                 transform=im2_ax.transAxes, ha='left', va='center', fontsize=11)
    fig.subplots_adjust(
    left=0.1, 
    right=0.85, 
    bottom=0.1,
    top=0.95, 
    wspace=0.15,
    hspace=0.4
    )
    fig.text(0.03, 0.7, r"$\mathrm{Counts\ (\times10^{6})}$", 
             ha='center', va='center', rotation='vertical', fontsize=16)
    #plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/spectra_360ks_logE.png")
    plt.show()
        
        
collapse_all(mass_output_array, output_exp_time)

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def collapse_all_separate(mass_output_array, output_exp_time,
                          zoom_center=20800, zoom_halfwidth=50,
                          pixel_scale_mas=7,
                          outdir="/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/360ks_spectra_logE"):
    """
    Makes one (spectrum + collapsed image) figure per mass and saves each individually.
    """

    os.makedirs(outdir, exist_ok=True)

    mpl.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 18
    })

    # If you have a 1-1 mapping between rows and these constants, keep it here.
    # Otherwise, better to compute log10(M) from mass_output_array directly.
    constants = np.array([ 6.0, 5.5, 5.0, 4.5, 4.0])

    for i, mass in enumerate(mass_output_array):
        base = f"{dir_basic}/{output_exp_time}ks_exposures/M_{mass}_output"

        data = np.load(f"{base}/data.npy")
        spectrum = np.load(f"{base}/spectrum_counts.npy")
        flux = np.load(f"{base}/spectrum_flux.npy")
        spectrum_std = np.load(f"{base}/spectrum_std.npy")
        wavelength_angstrom = np.load(f"{base}/wavelength_angstrom.npy")

        # --- Collapse cube around line ---
        mask = (wavelength_angstrom > zoom_center - zoom_halfwidth) & \
               (wavelength_angstrom < zoom_center + zoom_halfwidth)

        # NOTE: assumes data is (nwav, ny, nx)
        collapsed_image = np.nan_to_num(np.nansum(data[mask, :, :], axis=0))

        spectrum_cut = spectrum[mask]
        wavelength_cut = wavelength_angstrom[mask]

        # --- Make a single figure for this mass ---
        fig, axes = plt.subplots(
            1, 2, figsize=(10, 4),
            gridspec_kw={'width_ratios': [1.5, 1]}
        )

        #fig.subplots_adjust(wspace=0.35)
        spec_ax, im_ax = axes

        # Spectrum
        spec_ax.plot(wavelength_cut, spectrum_cut * 1e-6)
        try:
            popt, pcov = curve_fit(gaussian, wavelength_cut, spectrum_cut, p0=[1])
            spec_ax.plot(wavelength_cut, gaussian(wavelength_cut, *popt) * 1e-6,
                         ls='--', c='red')
        except Exception:
            # If the fit fails at low SNR, still save the plot without the model
            pass

        spec_ax.set_xlim(zoom_center - zoom_halfwidth, zoom_center + zoom_halfwidth)
        spec_ax.set_xlabel(r"$\lambda\ (\AA)$")
        spec_ax.set_ylabel(r"$\mathrm{Counts\ (}\times 10^{6}\mathrm{)}$")

        # Image
        ny, nx = collapsed_image.shape
        extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
        im = im_ax.imshow(collapsed_image, origin='lower', cmap='gray',
                          aspect='equal', extent=extent)
        im_ax.set_xlabel(r"$X\ (mas)$")
        im_ax.set_ylabel(r"$Y\ (mas)$")
        im_ax.yaxis.set_label_position("right")
        im_ax.yaxis.tick_right()
        im_ax.tick_params(labelleft=False)

        # Text annotation (mass + SNR)
        SNR_val, total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)

        # If your constants list matches your mass list, use it; otherwise compute from mass.
        if i < len(constants):
            logM = constants[i]
        else:
            # fallback: try to interpret `mass` as log10 mass if it's numeric
            try:
                logM = float(mass)
            except Exception:
                logM = np.nan

        fig.suptitle(
    fr"$M_{{\mathrm{{cluster}}}} = 10^{{{logM:.2f}}}\,M_{{\odot}}"
    r"\ \ \ \ \ "
    fr"\mathrm{{SNR}} = {SNR_val:.2f}$",
    fontsize=16, y=0.98
)

        fig.subplots_adjust(left=0.10, right=0.78, bottom=0.18, top=0.92, wspace=0.15)
        fig.subplots_adjust(wspace=0.35)
        fig.subplots_adjust(top=0.88)


        # Save
        outfile = os.path.join(outdir, f"spectra_{output_exp_time}ks_M_{mass}.png")
        #fig.savefig(outfile, dpi=200, bbox_inches="tight")
        plt.close(fig)

    print(f"Saved {len(mass_output_array)} figures to: {outdir}")


collapse_all_separate(mass_output_array, output_exp_time,
                          zoom_center=20800, zoom_halfwidth=50,
                          pixel_scale_mas=7,
                          outdir="/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/360ks_spectra_logE")

def display_all(mass_output_array, output_exp_time_array):
    fig, axes = plt.subplots( len(mass_output_array),len(output_exp_time_array) , figsize=(10, 35),
                              gridspec_kw={'width_ratios': [1, 1,1]}, sharey='row')
    #axes = axes.flatten()
    # row_ymins = [-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,-0.01,-0.005]
    # row_ymaxs = [0.2, 0.08, 0.02,0.008, 0.003, 0.002,0.002,0.002, 0.002]
    for t in range(len(output_exp_time_array)):
        for i in range(len(mass_output_array)):
            data = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/data.npy")
            spectrum = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/spectrum_counts.npy")
            flux = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/spectrum_flux.npy")
            spectrum_std = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/spectrum_std.npy")
            data_std = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/data_std.npy")
            data_flux = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/data_flux.npy")
            wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time_array[t]}ks_exposures/M_{mass_output_array[i]}_output/wavelength_angstrom.npy")
        
            zoom_center = 20800
            zoom_halfwidth = 50
            mask = (wavelength_angstrom > zoom_center - zoom_halfwidth) & (wavelength_angstrom < zoom_center + zoom_halfwidth)
            collapsed_image = np.nansum(data[mask,:,:],axis=0)
            collapsed_image = np.nan_to_num(collapsed_image)
            spectrum_cut = spectrum[mask]
            wavelength_cut = wavelength_angstrom[mask]

            n_rows, n_cols = axes.shape
            zoom  = axes[i,t]
            zoom.plot(wavelength_cut, spectrum_cut*10**(-6))
            popt, pcov = curve_fit(gaussian, wavelength_cut,spectrum_cut, p0=[1])
            zoom.plot(wavelength_cut, gaussian(wavelength_cut, *popt)*10**(-6), ls='--', c='red')
            zoom.set_xlim(zoom_center - zoom_halfwidth, zoom_center + zoom_halfwidth)
            #zoom.set_ylim(row_ymins[i], row_ymaxs[i])

            if i != n_rows - 1:
                zoom.set_xticklabels([])
            else:
                zoom.tick_params(axis='x',labelsize=8 ,labelbottom=True)
                plt.setp(zoom.get_xticklabels(), rotation=35, ha='right')

            
            # if t != 0:
            #     zoom.set_yticklabels([])
            # else:
            #     zoom.tick_params(axis='y', labelsize=8, labelleft=True)
            #     zoom.yaxis.set_major_locator(MaxNLocator(3))
            if t == 0:
                zoom.tick_params(axis='y', labelsize=8)
                zoom.yaxis.set_major_locator(MaxNLocator(3))

    
            if i == 0:
                zoom.set_title(f"$t_{{exp}} = {output_exp_time_array[t]} ks$", fontsize=16, pad=10)

            if t == n_cols - 1:
                zoom.text(1.05, 0.5,
                        rf"$10^{{{mass_output_array[i]}}}\ M_{{\odot}}$",
                        transform=zoom.transAxes,
                        va="center", ha="left",
                        fontsize=18)

    fig.subplots_adjust(
        left=0.08,
        right=0.86,
        bottom=0.12,
        top=0.93,
        wspace=0.15,
        hspace=0.25
    )

    fig.text(0.02, 0.7, r"Counts ($\times 10^{6}$)",
             ha='center', va='center', rotation="vertical", fontsize=20)

    fig.text(0.95, 0.05, r"$\lambda\ (\mathrm{\AA})$",
             ha='center', va='center', fontsize=20)
    #plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/data_logE.png")
    plt.show()
    

mass_output_array = [6.0,5.5,5.0,4.5]
output_exp_time_array =[360]
#collapse_all_separate(mass_output_array, output_exp_time,
#                          zoom_center=20800, zoom_halfwidth=50,
#                          pixel_scale_mas=7,
#                          outdir="/home/steff/hsim/zackrisson_pop3_all/code/Rep#ort_plots/interim_report/360ks_spectra_salpeter")


#collapse_all(mass_output_array, output_exp_time_array)






imf = ['salpeter','logA', 'logE']
output_exp_time =[360]
output_masses = [4.5,5.0,5.5,6.0]
shapes = ["^", "s", "o"]

ax = plt.subplot(111)
ax.set_yscale('log')
colors = cm.jet(np.linspace(0, 1, len(output_exp_time)))
for t in range(len(output_exp_time)):
    line_flux_array = []
    line_flux_err = []
    for i in range(len(output_masses)):
        spectrum = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_counts.npy")
        spectrum_flux = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_flux.npy")
        spectrum_std = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_std.npy")
        wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/wavelength_angstrom.npy")

        SNR, total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)
        line_flux, flux_err = calc_line_flux(spectrum_flux, spectrum_std, wavelength_angstrom)
        line_flux_array.append(line_flux)
        line_flux_err.append(flux_err)
        
    ax.scatter(output_masses, line_flux_array, color=colors[t],marker = shapes[t],s=40, label = f'\(t_{{exp}}: {output_exp_time[t]} ks\)')
    ax.errorbar(output_masses, line_flux_array, yerr=line_flux_err, fmt='none', color=colors[t])


plt.xlabel("\(M_{cluster}\ (10^{x \cdot M_{\odot}})\)")
plt.ylabel("\(F_{He_{II}\lambda 1640}\ (erg \cdot\ s^{-1}\cdot cm^{-2})\)")
plt.legend()
#plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/mass_flux.png")
plt.show()
output_hours = [3, 14, 100]
output_exp_time = [10,50,360]
ax = plt.subplot(111)
plt.yscale('log', base=10)
plt.xscale('log', base=10)
colors = cm.jet(np.linspace(0, 1, len(output_exp_time)))
for t in range(len(output_exp_time)):
    flux_array = []
    flux_error = []
    SNR_array = []
    for i in range(len(output_masses)):
        spectrum = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_counts.npy")
        spectrum_flux = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_flux.npy")
        spectrum_std = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_std.npy")
        wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/wavelength_angstrom.npy")

        SNR, total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)
        line_flux, flux_err = calc_line_flux(spectrum_flux, spectrum_std, wavelength_angstrom)
        flux_array.append(line_flux)
        flux_error.append(flux_err)
        SNR_array.append(SNR)
    plt.scatter(flux_array, SNR_array,color = colors[t], marker = shapes[t], s=40,label = f'\(t_{{exp}}: {output_exp_time[t]} ks\ (\sim {output_hours[t]} h)\)')
    plt.errorbar(flux_array, SNR_array, xerr=line_flux_err, fmt='none', color=colors[t])


plt.ylabel("\(SNR\)")
plt.xlabel("\(F_{He_{II}\lambda 1640}\ (erg \cdot\ s^{-1}\cdot cm^{-2})\)")
plt.legend()
#plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/SNR_flux_logA.png")
plt.show()


output_exp_time =[360]

output_masses = np.array([4.5,5.0,5.5,6.0])
output_masses_10s = 10**(output_masses)

#plt.figure(figsize=(8, 8))
ax = plt.subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
colors = ['#0072B2','#E69F00','#009E73']
shapes = ["^", "s", "o"]
symmetric_perc = []
for t in range(len(output_exp_time)):
    line_flux_array = []
    input_line_flux_array = []
    line_flux_array_err = []
    SNR_array = []
    SNR_2_array = []
    symmetric_perc.append([])
    for i in range(len(output_masses)):
        spectrum = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_counts.npy")
        spectrum = spectrum / 4*np.pi
        spectrum_flux = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_flux.npy")
        spectrum_flux = spectrum_flux / 4*np.pi
        spectrum_std = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_std.npy")
        wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/wavelength_angstrom.npy")

        input_flux = np.load(f"/home/steff/hsim/HSIM/hsim/input_cubes/salpeter/360ks_exposures/M_{output_masses[i]}_flux.npy")
        input_wavelength = np.load(f"/home/steff/hsim/HSIM/hsim/input_cubes/salpeter/360ks_exposures/M_{output_masses[i]}_wavelength.npy")
        input_std = np.zeros(len(input_wavelength))

        SNR, total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)
        SNR_2 = calcSNR_peak(wavelength_angstrom, spectrum, spectrum_std)
        line_flux, flux_err = calc_line_flux(spectrum_flux, spectrum_std, wavelength_angstrom)
        input_line_flux, input_err = calc_line_flux(input_flux, input_std, input_wavelength)
        line_flux_array.append(line_flux)
        input_line_flux_array.append(input_line_flux)
        line_flux_array_err.append(flux_err)
        SNR_array.append(SNR)
        SNR_2_array.append(SNR_2)
        delta_SNR = ((SNR_2 - SNR)/SNR) * 100
        symmetric_perc[t].append(delta_SNR)
        
    print(f"Methods differ by {symmetric_perc}%")
    ax.scatter(output_masses_10s, SNR_array,marker = shapes[t] ,color=colors[t],s=30, label = fr'\(t_{{exp}}\ = {output_exp_time[t]}\ ks\)')# (\approx {output_hours[t]} h)\)')
    ax.plot(output_masses_10s, SNR_array ,color=colors[t])
    #ax.scatter(output_masses, SNR_2_array, marker = shapes[t],s=70, color = colors[t],edgecolors='black')

ax.hlines(y = 4.9, xmin=10**(3.5), xmax=7*10**(4), ls='--', color='#009E73')
ax.hlines(y = 5, xmin=10**(3.5), xmax=2*10**(5), ls='--', color='#E69F00')
ax.hlines(y = 5.1, xmin=10**(3.5), xmax=3*10**(5), ls='--', color='#0072B2')
ax.vlines(x = 7*10**4, ymin = 0, ymax = 5, ls = '--', color = '#009E73')
ax.vlines(x = 2*10**5, ymin = 0, ymax = 5, ls = '--', color = '#E69F00')
ax.vlines(x = 3*10**5, ymin = 0, ymax = 5, ls = '--', color = '#0072B2')
plt.xlim(10**(3.5), 10**(8.5))
plt.xlabel("\(M_{cluster}\ (M_{\odot})\)")
plt.ylabel("\(SNR\)")
plt.text(1*10**4, 4*10**2, '$IMF: Salpeter$\n$z = 11.5$', bbox=dict(
        boxstyle='round',
        facecolor='white',
        alpha=0.8
    ))
plt.legend(loc = 'upper center',ncols = 3,  bbox_to_anchor =(0.5, 1.19),frameon = False,
    borderpad=0.2,
    labelspacing=0.3,
    handletextpad=0.4
)

#ax.add_artist(legend_main)
# no1 = ax.scatter([], [],                   
#            marker='o',
#            edgecolors='black',
#            facecolors='none',
#            linewidths=1.2,
#            label='\(Gaussian\ peak\ method\)')
# no2 = ax.scatter([],[],
#            marker = 'o', label='\(\chi^2 method\)')
# plt.legend(handles=[no1,no2], loc = 'lower right')
plt.xlim(8*10**(3), 10**(8))
#plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/SNR_mass_sal.png")
plt.show()


# --- NEW PLOT: test SNR ∝ sqrt(t) by checking SNR/sqrt(t) ---
fig, ax = plt.subplots()

ax.set_xscale('log')

colors = ['#0072B2','#E69F00','#009E73']
shapes = ["^", "s", "o"]

for t_idx, t_ks in enumerate(output_exp_time):
    SNR_array = []
    
    for i in range(len(output_masses)):
        spectrum = np.load(f"{dir_basic}/{t_ks}ks_exposures/M_{output_masses[i]}_output/spectrum_counts.npy")
        spectrum_std = np.load(f"{dir_basic}/{t_ks}ks_exposures/M_{output_masses[i]}_output/spectrum_std.npy")
        wavelength_angstrom = np.load(f"{dir_basic}/{t_ks}ks_exposures/M_{output_masses[i]}_output/wavelength_angstrom.npy")

        SNR, _ = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)
        SNR_array.append(SNR)

    # IMPORTANT: use seconds (or ks) consistently; sqrt scaling works with any time unit
    snr_scaled = np.array(SNR_array) / np.sqrt(t_ks)  # using ks units is fine for a proportionality check

    ax.scatter(output_masses_10s, snr_scaled, marker=shapes[t_idx], color=colors[t_idx],
            label=fr"$t_{{exp}}={t_ks}\ \mathrm{{ks}}$")

    ax.plot(output_masses_10s, snr_scaled, marker='none', color=colors[t_idx])

ax.set_xlabel(r"$M_{\rm cluster}\ (M_\odot)$")
ax.set_ylabel(r"$\mathrm{SNR}/\sqrt{t_{\rm exp}\ (\mathrm{ks})}$")
ax.legend(loc = 'upper center',ncols = 3,  bbox_to_anchor =(0.5, 1.15),frameon = False,
    borderpad=0.2,
    labelspacing=0.3,
    handletextpad=0.4
)

#plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/SNR_sqrt_t_check_logA.png")
plt.show()

SNR_target = 5.0
t0 = 360  # ks baseline

# Compute SNR at t0 for each mass
SNR_t0_sal = []

for m in output_masses:
    spectrum = np.load(f"{dir_basic}/{t0}ks_exposures/M_{m}_output/spectrum_counts.npy")
    spectrum_std = np.load(f"{dir_basic}/{t0}ks_exposures/M_{m}_output/spectrum_std.npy")
    wave = np.load(f"{dir_basic}/{t0}ks_exposures/M_{m}_output/wavelength_angstrom.npy")
    snr, _ = calcSNR_chi2(wave, spectrum, spectrum_std)
    SNR_t0_sal.append(snr)


SNR_t0_sal = np.array(SNR_t0_sal)

mask = SNR_t0_sal > 0.5
masses_lin = 10**output_masses[mask]
t_req_sal = t0 * (SNR_target / SNR_t0_sal[mask])**2  # ks
M_ref = np.logspace(4, 6.3, 200)
M0 = 1*10**(5.5)
texp = t_req_sal[3]
t_ref = texp * (M_ref / M0)**(-2)


plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.plot(masses_lin, t_req_sal, marker='o', ms = 5, lw = 0.5)
plt.plot(
    M_ref, t_ref,
    linestyle='--',
    color='red',
    linewidth=1,
    label=r'$t_{\rm exp} \propto M_{\rm cluster}^{-2}$'
)
plt.xlabel(r"$M_{\rm cluster}\ (M_\odot)$")
plt.ylabel(r"$min^{m}\ t_{\rm exp}\ (\mathrm{ks})$")
plt.xlim(8*10**3,2*10**6 )
plt.ylim(10**(-1), 10**4)
plt.legend(loc = 'upper center',ncols = 2,  bbox_to_anchor =(0.5, 1.19),frameon = False,
    borderpad=0.2,
    labelspacing=0.3,
    handletextpad=0.4
)
plt.text(1*10**4, 0.5, '$IMF: Salpeter$\n$z = 11.5$', bbox=dict(
        boxstyle='round',
        facecolor='white',
        alpha=0.8
    ))
#plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/SNR_texp_logA.png")
plt.show()



# Compute SNR at t0 for each mass
SNR_t0_sal = []
SNR_t0_logA = []
SNR_t0_logE = []

for m in output_masses:
    spectrum = np.load(f"{dir_basic}/{t0}ks_exposures/M_{m}_output/spectrum_counts.npy")
    spectrum = spectrum / 4*np.pi
    print(f"spectrum salpeter: {spectrum}")
    spectrum_std = np.load(f"{dir_basic}/{t0}ks_exposures/M_{m}_output/spectrum_std.npy")
    wave = np.load(f"{dir_basic}/{t0}ks_exposures/M_{m}_output/wavelength_angstrom.npy")
    snr, _ = calcSNR_chi2(wave, spectrum, spectrum_std)
    SNR_t0_sal.append(snr)

for m in output_masses:
    spectrum = np.load(f"{dir_basic_logA}/{t0}ks_exposures/M_{m}_output/spectrum_counts.npy")
    spectrum = spectrum / 4*np.pi
    spectrum_std = np.load(f"{dir_basic_logA}/{t0}ks_exposures/M_{m}_output/spectrum_std.npy")
    wave = np.load(f"{dir_basic_logA}/{t0}ks_exposures/M_{m}_output/wavelength_angstrom.npy")
    snr, _ = calcSNR_chi2(wave, spectrum, spectrum_std)
    SNR_t0_logA.append(snr)

for m in output_masses:
    spectrum = np.load(f"{dir_basic_logE}/{t0}ks_exposures/M_{m}_output/spectrum_counts.npy")
    spectrum = spectrum / 4*np.pi
    print(f"spectrum logE: {spectrum}")
    spectrum_std = np.load(f"{dir_basic_logE}/{t0}ks_exposures/M_{m}_output/spectrum_std.npy")
    wave = np.load(f"{dir_basic_logE}/{t0}ks_exposures/M_{m}_output/wavelength_angstrom.npy")
    snr, _ = calcSNR_chi2(wave, spectrum, spectrum_std)
    SNR_t0_logE.append(snr)

SNR_t0_sal = np.array(SNR_t0_sal)
SNR_t0_logA = np.array(SNR_t0_logA)
SNR_t0_logE = np.array(SNR_t0_logE)

mask = SNR_t0_sal > 0.5
masses_lin = 10**output_masses[mask]
t_req_sal = t0 * (SNR_target / SNR_t0_sal[mask])**2  # ks
M_ref = np.logspace(3, 6.3, 200)
M0 = 1*10**(5.5)
texp = t_req_sal[2]
t_ref_sal = texp * (M_ref / M0)**(-2)

mask = SNR_t0_logA > 0.5
masses_lin_logA = 10**output_masses[mask]
t_req_logA = t0 * (SNR_target / SNR_t0_logA[mask])**2  # ks
M_ref = np.logspace(3, 6.3, 200)
M0 = 1*10**(5.5)
texp = t_req_logA[2]
t_ref_logA = texp * (M_ref / M0)**(-2)

mask = SNR_t0_logE > 0.5
masses_lin_logE = 10**output_masses[mask]
t_req_logE = t0 * (SNR_target / SNR_t0_logE[mask])**2  # ks
M_ref = np.logspace(3, 6.3, 200)
M0 = 1*10**(5.5)
texp = t_req_logE[2]
t_ref_logE = texp * (M_ref / M0)**(-2)

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.scatter(masses_lin, t_req_sal, marker='o', label = '$Salpeter$')
plt.scatter(masses_lin_logA, t_req_logA, marker='o', label = '$logA$')
plt.scatter(masses_lin_logE, t_req_logE, marker='o', label = '$logE$')
plt.plot(
    M_ref, t_ref_sal,
    linestyle='--',
    color='blue',
    linewidth=1,
    alpha = 0.6
    )
plt.plot(
    M_ref, t_ref_logA,
    linestyle='--',
    color='green',
    linewidth=1,
    alpha = 0.6
    )
plt.plot(
    M_ref, t_ref_logE,
    linestyle='--',
    color='orange',
    linewidth=1,
    alpha = 0.6
    )
plt.xlabel(r"$M_{\rm cluster}\ (M_\odot)$")
plt.ylabel(r"$min^{m}\ t_{\rm exp}\ (\mathrm{ks})$")
plt.xlim(5*10**3,5*10**5 )
plt.ylim(3*10**(0), 3*10**3)
plt.legend(loc = 'upper center',ncols = 3,  bbox_to_anchor =(0.5, 1.19),frameon = False,
    borderpad=0.2,
    labelspacing=0.3,
    handletextpad=0.4
)
plt.text(9*10**4,10**3, '$z = 11.5$', bbox=dict(
        boxstyle='round',
        facecolor='white',
        alpha=0.8
    ))
plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/SNR_texp_logE.png")
plt.show()






ax = plt.subplot(111)
ax.scatter(output_masses[3:],symmetric_perc[0][3:], marker = shapes[0], color = colors[0], label = f'$t_{{exp}} = {output_exp_time[0]}s$')
ax.scatter(output_masses[2:],symmetric_perc[1][2:], marker = shapes[1], color = colors[1], label = f'$t_{{exp}} = {output_exp_time[1]}s$')
ax.scatter(output_masses[0:],symmetric_perc[2][0:], marker = shapes[2], color = colors[2], label = f'$t_{{exp}} = {output_exp_time[2]}s$')
ax.legend()
plt.show()

ax = plt.subplot(111)
ax.set_yscale('log')
colors = cm.jet(np.linspace(0, 1, len(output_exp_time)))
shapes = ["^", "s", "o"]
for t in range(len(output_exp_time)):
    line_flux_array = []
    line_flux_array_err = []
    input_line_flux_array = []
    SNR_array = []
    for i in range(len(output_masses)):
        spectrum = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_counts.npy")
        spectrum_flux = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_flux.npy")
        spectrum_std = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/spectrum_std.npy")
        wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time[t]}ks_exposures/M_{output_masses[i]}_output/wavelength_angstrom.npy")
        input_flux = np.load(f"/home/steff/hsim/HSIM/hsim/input_cubes/salpeter/10ks_exposures/M_{output_masses[i]}_flux.npy")
        input_wavelength = np.load(f"/home/steff/hsim/HSIM/hsim/input_cubes/salpeter/10ks_exposures/M_{output_masses[i]}_wavelength.npy")
        input_std = np.zeros(len(input_wavelength))

        SNR, total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_ge0_sal_500_001_is5.22"
        data = ascii.read(file_loc,guess = True, data_start = 0)
        H_beta = data['col5']
        HeII_1640 = data['col15'] * H_beta
        z = 11.5
        comoving_d = cosmo.comoving_distance(z)
        d_lumo_cm = (1 + z)* (comoving_d * (3.086*10**24))
        input_line_flux_prez = (HeII_1640[0])#/(np.sqrt(2*np.pi)*sigma_line))
        input_line_flux = input_line_flux_prez * (1/(d_lumo_cm**2)) * (1/(1 + z))
        input_line_flux = input_line_flux.value
        input_line_flux = input_line_flux * 10**(output_masses[i])
        input_line_flux = input_line_flux * (1+z)#* (input_scale)**2
        line_flux, flux_err = calc_line_flux(spectrum_flux, spectrum_std, wavelength_angstrom)
        #input_line_flux, input_err = calc_line_flux(input_flux, input_std, input_wavelength)
        line_flux_array.append(line_flux)
        input_line_flux_array.append(input_line_flux)
        line_flux_array_err.append(flux_err)
        SNR_array.append(SNR)
    print(input_line_flux_array)
    print(line_flux_array)
    ratio = np.array(input_line_flux_array) / np.array(line_flux_array)
    error = ratio * np.sqrt(0 + (flux_err / line_flux)**2)
    error = np.abs(error)
    ax.scatter(output_masses, ratio, color=colors[t],marker = shapes[t] ,s=40,label = f'\(t_{{exp}}: {output_exp_time[t]} ks\)')
    ax.errorbar(output_masses, ratio, yerr=error, fmt='none', color=colors[t])

plt.xlabel("\(M_{cluster}\ (10^{x \cdot M_{\odot}})\)")
plt.ylabel("\(f_{input} / f_{output}\)")
plt.legend()
#plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/mass_ratio_log_logA.png")
plt.show()



