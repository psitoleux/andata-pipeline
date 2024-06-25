from abmodule import Module
import scanpy as sc
import os
import subprocess




def log1p(adata : sc.AnnData, L : int = 1e4) -> sc.AnnData:
    sc.pp.normalize_total(adata, target_sum=L)
    sc.pp.log1p(adata)
    
    return adata


# TODO GPT-4o generated function, clean it up and understand it !!!!

def sanity_normalization(
    adata: sc.AnnData,
    sanity_executable: str = "path_to_Sanity_binary",
    destination_folder: str = ".",
    n_threads: int = 4,
    extended_output: bool = False,
    vmin: float = 0.001,
    vmax: float = 50,
    nbin: int = 160,
    no_cell_size_normalization: bool = False
) -> sc.AnnData:
    """
    Apply Sanity normalization to an AnnData object.

    Parameters:
    - adata: AnnData object with raw UMI counts in adata.X.
    - sanity_executable: Path to the Sanity executable.
    - destination_folder: Destination folder for Sanity output files.
    - n_threads: Number of threads to use.
    - extended_output: Whether to generate extended output.
    - vmin: Minimal value of variance in log transcription quotient.
    - vmax: Maximal value of variance in log transcription quotient.
    - nbin: Number of bins for the variance in log transcription quotient.
    - no_cell_size_normalization: Whether to skip cell size normalization.

    Returns:
    - AnnData object with normalized data.
    """
    # Ensure the destination folder exists
    os.makedirs(destination_folder, exist_ok=True)

    # Save the raw UMI count matrix to a temporary file
    umi_count_file = os.path.join(destination_folder, "umi_counts.txt")
    np.savetxt(umi_count_file, adata.X.T, delimiter='\t')

    # Prepare the Sanity command
    sanity_cmd = [
        sanity_executable,
        "-f", umi_count_file,
        "-d", destination_folder,
        "-n", str(n_threads),
        "-vmin", str(vmin),
        "-vmax", str(vmax),
        "-nbin", str(nbin)
    ]

    if extended_output:
        sanity_cmd.extend(["-e", "1"])
    if no_cell_size_normalization:
        sanity_cmd.extend(["-no_norm", "1"])

    # Run Sanity normalization
    subprocess.run(sanity_cmd, check=True)

    # Load the normalized data
    normalized_data_file = os.path.join(destination_folder, "log_transcription_quotients.txt")
    normalized_data = np.loadtxt(normalized_data_file, delimiter='\t', skiprows=1, usecols=range(1, adata.shape[0] + 1)).T

    # Update the AnnData object with normalized data
    adata.X = normalized_data

    # Optionally, load extended output if required
    if extended_output:
        mu_file = os.path.join(destination_folder, "mu.txt")
        variance_file = os.path.join(destination_folder, "variance.txt")
        delta_file = os.path.join(destination_folder, "delta.txt")

        adata.var["mu"] = np.loadtxt(mu_file, delimiter='\t', skiprows=1, usecols=1)
        adata.var["variance"] = np.loadtxt(variance_file, delimiter='\t', skiprows=1, usecols=1)
        adata.layers["delta"] = np.loadtxt(delta_file, delimiter='\t', skiprows=1, usecols=range(1, adata.shape[0] + 1)).T

    return adata