import ROOT
import numpy as np

def process_histogram(input_file, hist_name, output_file):
    # Open the input ROOT file
    file = ROOT.TFile(input_file, "READ")
    hist = file.Get(hist_name)

    if not hist:
        print(f"Histogram {hist_name} not found in file {input_file}")
        return

    # Clone the histogram to modify it
    hist_clone = hist.Clone(f"{hist.GetName()}")

    # Get the number of bins in x and y directions
    nx_bins = hist.GetNbinsX()
    ny_bins = hist.GetNbinsY()

    # Loop over all bins
    for x in range(1, nx_bins + 1):
        for y in range(1, ny_bins + 1):
            bin_content = hist.GetBinContent(x, y)
            if bin_content < 0:
                # Calculate the average of surrounding bins
                surrounding_bins = []
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        if dx == 0 and dy == 0:
                            continue
                        nx, ny = x + dx, y + dy
                        if 1 <= nx <= nx_bins and 1 <= ny <= ny_bins:
                            surrounding_bins.append(hist.GetBinContent(nx, ny))

                # Replace the bin content with the average of the surrounding bins
                if surrounding_bins:
                    avg_value = np.mean(surrounding_bins)
                else:
                    avg_value = 0  # Default value if no surrounding bins
                hist_clone.SetBinContent(x, y, avg_value)

    # Save the modified histogram to a new ROOT file
    output = ROOT.TFile(output_file, "RECREATE")
    hist_clone.Write()
    output.Close()
    file.Close()

    print(f"Processed histogram saved to {output_file}")

# Example usage
process_histogram("/vols/cms/ks1021/TIDAL/Draw/plots/zpt_reweighting_LO_2023BPix.root", "zptmass_histo", "/vols/cms/ks1021/TIDAL/Draw/plots/zpt_reweighting_LO_2023BPix_no_negative_bins.root")
