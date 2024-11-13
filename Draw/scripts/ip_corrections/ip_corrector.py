import ROOT
import numpy as np

def load_histograms(file_name, channel, ip_names, ip_err_names, ip_sig_names):

    file = ROOT.TFile.Open(file_name)
    if not file or file.IsZombie():
        raise FileNotFoundError(f"IpCorrection: file {file_name} does not exist or is corrupted")

    hist_names = [key.GetName() for key in file.GetListOfKeys()]

    eta_names = list(set([name.split("_")[-2] for name in hist_names if name.endswith("_data")]))
    eta_names = sorted(eta_names)

    eta_bin_edges = sorted([name.replace("p", ".").replace("to", "-") for name in eta_names])
    eta_ranges = [float(i.split("-")[0]) for i in eta_bin_edges]
    eta_ranges.append(float(eta_bin_edges[-1].split("-")[1]))

    # Initialize empty dictionaries to store histograms by category
    hist_ip_data = {ip_name: {eta_name: None for eta_name in eta_names} for ip_name in ip_names}
    hist_ip_mc = {ip_name: {eta_name: None for eta_name in eta_names} for ip_name in ip_names}
    hist_err_data = {ip_err_name: {eta_name: None for eta_name in eta_names} for ip_err_name in ip_err_names}
    hist_err_mc = {ip_err_name: {eta_name: None for eta_name in eta_names} for ip_err_name in ip_err_names}
    hist_sig_data = {ip_sig_name: {eta_name: None for eta_name in eta_names} for ip_sig_name in ip_sig_names}
    hist_sig_mc = {ip_sig_name: {eta_name: None for eta_name in eta_names} for ip_sig_name in ip_sig_names}

    # Helper function to convert ROOT histogram to (counts, edges) numpy arrays
    def convert_hist_to_numpy(hist):
        bin_edges = np.array([hist.GetBinLowEdge(i+1) for i in range(hist.GetNbinsX() + 1)])
        bin_counts = np.array([hist.GetBinContent(i+1) for i in range(hist.GetNbinsX())])
        return bin_counts, bin_edges

    # Helper function to load histograms and convert them to numpy arrays
    def load_hist_as_numpy(hist_name):
        hist = file.Get(hist_name)
        if hist:
            return convert_hist_to_numpy(hist)
        else:
            print(f"Warning: histogram {hist_name} does not exist")
        return None

    # Load and convert IP histograms
    for ip_name in ip_names:
        for eta_name in eta_names:
            hist_ip_data[ip_name][eta_name] = load_hist_as_numpy(f"{ip_name}_{eta_name}_data")
            hist_ip_mc[ip_name][eta_name] = load_hist_as_numpy(f"{ip_name}_{eta_name}_MC")

    # Load and convert IP error histograms
    for ip_err_name in ip_err_names:
        for eta_name in eta_names:
            hist_err_data[ip_err_name][eta_name] = load_hist_as_numpy(f"{ip_err_name}_{eta_name}_data")
            hist_err_mc[ip_err_name][eta_name] = load_hist_as_numpy(f"{ip_err_name}_{eta_name}_MC")

    # Load and convert IP significance histograms
    for ip_sig_name in ip_sig_names:
        for eta_name in eta_names:
            hist_sig_data[ip_sig_name][eta_name] = load_hist_as_numpy(f"{ip_sig_name}_{eta_name}_data")
            hist_sig_mc[ip_sig_name][eta_name] = load_hist_as_numpy(f"{ip_sig_name}_{eta_name}_MC")

    return {
        "eta_names": eta_names,
        "eta_ranges": eta_ranges,
        "hist_ip_data": hist_ip_data,
        "hist_ip_mc": hist_ip_mc,
        "hist_err_data": hist_err_data,
        "hist_err_mc": hist_err_mc,
        "hist_sig_data": hist_sig_data,
        "hist_sig_mc": hist_sig_mc
    }


def correct_IP(hist_data, ip_x, ip_y, ip_z, gen_ip_x, gen_ip_y, gen_ip_z, eta):
    # Determine the eta bin for each value in eta array based on absEta
    absEta = np.abs(eta)
    eta_bins = np.digitize(absEta, bins=hist_data["eta_ranges"]) - 1  # Find eta bin for each entry
    eta_bins = np.clip(eta_bins, 0, len(hist_data["eta_names"]) - 1)  # Ensure indices are within valid range

    # Helper function to apply quantile mapping for a single component
    def apply_mapping(ip, eta_bin_index, component):
        eta_name = hist_data["eta_names"][eta_bin_index]
        hist_MC = hist_data["hist_ip_mc"][component][eta_name]
        hist_data_ = hist_data["hist_ip_data"][component][eta_name]
        return applyQuantileMapping_numpy(hist_MC[0], hist_MC[1], hist_data_[0], hist_data_[1], ip)

    # Vectorized quantile mapping for each IP component
    apply_mapping_x = np.vectorize(apply_mapping, excluded=['component'])
    apply_mapping_y = np.vectorize(apply_mapping, excluded=['component'])
    apply_mapping_z = np.vectorize(apply_mapping, excluded=['component'])

    # Apply quantile mapping for each component based on the eta bin
    ip_x_corr = apply_mapping_x(ip_x-gen_ip_x, eta_bins, component="ip_x")
    ip_y_corr = apply_mapping_y(ip_y-gen_ip_y, eta_bins, component="ip_y")
    ip_z_corr = apply_mapping_z(ip_z-gen_ip_z, eta_bins, component="ip_z")

    ip_x_corr += gen_ip_x
    ip_y_corr += gen_ip_y
    ip_z_corr += gen_ip_z

    return ip_x_corr, ip_y_corr, ip_z_corr


def correct_ip_cov(hist_data, ip_cov_xx, ip_cov_yy, ip_cov_zz, eta):
    # Compute absolute value of eta and determine the appropriate bin indices for each entry
    abs_eta = np.abs(eta)
    eta_bins = np.digitize(abs_eta, bins=hist_data["eta_ranges"]) - 1
    eta_bins = np.clip(eta_bins, 0, len(hist_data["eta_names"]) - 1)  # Ensure indices are within bounds

    # Helper function to compute the corrected covariance for a single entry
    def compute_corrected_cov(eta_bin_idx, cov_xx, cov_yy, cov_zz):
        # Retrieve the correct eta name for the current bin
        eta_name = hist_data["eta_names"][eta_bin_idx]

        # Corrected error array for each component
        err = np.zeros(3)
        ip_cov_diag = [cov_xx, cov_yy, cov_zz]
        components = ["ip_x_Err", "ip_y_Err", "ip_z_Err"]

        for i, (component, cov_ii) in enumerate(zip(components, ip_cov_diag)):
            # Access the corresponding MC and data histograms for the component and eta bin
            hist_mc = hist_data["hist_err_mc"][component][eta_name]
            hist_data_ = hist_data["hist_err_data"][component][eta_name]

            # Calculate the standard deviation (sqrt of covariance) for the current component
            std_dev = np.sqrt(cov_ii)

            # Apply quantile mapping to map std_dev from MC to data histogram
            err[i] = applyQuantileMapping_numpy(hist_mc[0], hist_mc[1], hist_data_[0], hist_data_[1], std_dev)

            # Normalize the corrected error by the original standard deviation
            err[i] /= std_dev

        # Construct the corrected covariance matrix using the computed errors
        ip_cov_corrected = np.array([
            [err[0] * err[0] * cov_xx, err[0] * err[1] * np.sqrt(cov_xx * cov_yy), err[0] * err[2] * np.sqrt(cov_xx * cov_zz)],
            [err[1] * err[0] * np.sqrt(cov_yy * cov_xx), err[1] * err[1] * cov_yy, err[1] * err[2] * np.sqrt(cov_yy * cov_zz)],
            [err[2] * err[0] * np.sqrt(cov_zz * cov_xx), err[2] * err[1] * np.sqrt(cov_zz * cov_yy), err[2] * err[2] * cov_zz]
        ])

        return ip_cov_corrected

    # Vectorize compute_corrected_cov function to apply it element-wise across the input arrays
    vectorized_compute_corrected_cov = np.vectorize(compute_corrected_cov, signature='(),(),(),()->(3,3)')

    # Apply the vectorized function to each entry in eta
    ip_cov_corrected_all = vectorized_compute_corrected_cov(eta_bins, ip_cov_xx, ip_cov_yy, ip_cov_zz)

    return ip_cov_corrected_all


def applyQuantileMapping_numpy(hist_MC_counts, hist_MC_edges, hist_data_counts, hist_data_edges, var):
    # Define the range for `var` within the MC histogram
    xmin = hist_MC_edges[1]
    xmax = hist_MC_edges[-1]

    # Normalize both histograms
    normMC = np.sum(hist_MC_counts)
    normData = np.sum(hist_data_counts)

    # If `var` is outside the bounds of `hist_MC`, return `var` as is
    if var < xmin or var > xmax:
        return var

    # Find the bin in hist_MC that contains `var`
    nBinMC = np.searchsorted(hist_MC_edges, var, side="right") - 1

    # Compute cumulative probabilities at the lower and upper edges of the bin in hist_MC
    int_lower = np.sum(hist_MC_counts[:nBinMC]) / normMC
    int_upper = np.sum(hist_MC_counts[:nBinMC + 1]) / normMC
    int_diff = int_upper - int_lower

    # Get the bin width and lower edge in the MC histogram for interpolation
    xlow = hist_MC_edges[nBinMC]
    binWidth = hist_MC_edges[nBinMC + 1] - xlow

    # Compute the precise cumulative probability for `var` within the MC histogram bin
    int_center = int_lower + int_diff * (var - xlow) / binWidth

    # Find the corresponding bin in hist_data that matches the cumulative probability `int_center`
    int_data_lower = 0
    int_data_higher = 0
    nBinData = 0
    for j in range(len(hist_data_counts)):
        int_data_lower = np.sum(hist_data_counts[:j]) / normData
        int_data_higher = np.sum(hist_data_counts[:j + 1]) / normData
        if int_center > int_data_lower and int_center < int_data_higher:
            nBinData = j
            break

    # Get the bin width and lower edge in the data histogram for interpolation
    binWidthData = hist_data_edges[nBinData + 1] - hist_data_edges[nBinData]
    xlowData = hist_data_edges[nBinData]

    # Calculate the corrected value in the data histogram that corresponds to `int_center`
    int_data_diff = int_data_higher - int_data_lower
    varcorr = xlowData + binWidthData * (int_center - int_data_lower) / int_data_diff

    return varcorr


def corrected_IPSig(ipx, ipy, ipz, ip_cov):

    ip_vectors = np.vstack((ipx, ipy, ipz)).T
    ip_magnitudes = np.linalg.norm(ip_vectors, axis=1)
    ip_normalized = ip_vectors / ip_magnitudes[:, np.newaxis]

    error_ip = np.sqrt(np.einsum('ni,nij,nj->n', ip_normalized, ip_cov, ip_normalized))

    ip_significance = ip_magnitudes / error_ip

    return ip_significance

# Usage example:
path = "/vols/cms/ks1021/TIDAL/Draw/scripts/ip_corrections/ip_2018_pvbs.root"
path_2 = "/vols/cms/ks1021/TIDAL/Draw/plots/production_11_11_2024/Run3_2022/ip_calculation/histograms/combined_histograms_mm_Run3_2022.root"
file_name = path_2
ip_names = ["ip_x", "ip_y", "ip_z"]
ip_err_names = ["ip_x_Err", "ip_y_Err", "ip_z_Err"]
ip_sig_names = ["ip_LengthSig"]

# Check within function scope
hist_data = load_histograms(file_name, "mm", ip_names, ip_err_names, ip_sig_names)
# Example input arrays for testing
ip_x = np.array([0.001, 0.0015, 0.002, 0.0025])       # Example IP values for x
ip_y = np.array([0.002, 0.0025, 0.003, 0.0035])       # Example IP values for y
ip_z = np.array([0.003, 0.0035, 0.004, 0.0045])       # Example IP values for z

gen_ip_x = np.array([0.4, 0.45, 0.5, 0.55])   # Generated IP values for x
gen_ip_y = np.array([0.5, 0.55, 0.6, 0.65])   # Generated IP values for y
gen_ip_z = np.array([0.6, 0.65, 0.7, 0.75])   # Generated IP values for z

eta = np.array([0.91, 1.1, 1.9, 2.09])         # Eta values

ip_cov_xx = np.array([0.0001, 0.0002, 0.0003, 0.0004])  # Covariance values for IP x
ip_cov_yy = np.array([0.0002, 0.0003, 0.0004, 0.0005])  # Covariance values for IP y
ip_cov_zz = np.array([0.0003, 0.0004, 0.0005, 0.0006])  # Covariance values for IP z


# Run the correct_IP function
ip_x_corr, ip_y_corr, ip_z_corr = correct_IP(hist_data, ip_x, ip_y, ip_z, gen_ip_x, gen_ip_y, gen_ip_z, eta)

ip_cov = correct_ip_cov(hist_data, ip_cov_xx, ip_cov_yy, ip_cov_zz, eta)

ip_sig = corrected_IPSig(ip_x_corr, ip_y_corr, ip_z_corr, ip_cov)

# Print the corrected IP values
print("Corrected ip_x:", ip_x_corr)
print("Corrected ip_y:", ip_y_corr)
print("Corrected ip_z:", ip_z_corr)

# Print the corrected covariance matrices
print("Corrected covariance matrices:")
print(ip_cov)

print("Corrected IP significance:", ip_sig)
