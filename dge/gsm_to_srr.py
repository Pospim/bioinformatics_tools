"""
GSM to SRR Converter Script

This script retrieves SRA Run IDs (SRR) from the Gene Expression Omnibus (GEO)
database based on a list of GEO Sample IDs (GSM). It fetches the associated
SRX IDs and then retrieves the corresponding SRR IDs from the Sequence Read Archive
(SRA) database. The results are saved back into the original CSV file and also output list
as a separate text file for further use.

Usage:
    python gsm_to_srr.py /path/to/your/directory

Parameters:
    - workdir (str): The directory containing the input CSV file named 'samples.csv',
                     which should have a column labeled 'GSM' with the GSM IDs to be
                     processed.

Input File Format:
    The input CSV file (samples.csv) should have the following format:
    ```
    GSM
    GSM_ID_1
    GSM_ID_2
    ...
    ```

Output:
    - The original CSV file will be updated with a new column 'SRR', containing
      the corresponding SRR IDs for each GSM ID.
    - A new text file named 'SRRs.txt' will be created in the same directory,
      containing a bash-style list of all SRR IDs found.

Functions:
    - gsm_to_srx(gsm_id): Fetches the SRX ID associated with a given GSM ID.
    - srx_to_srr(srx_id): Fetches the SRR ID associated with a given SRX ID.
    - gsm_to_srx_batch(gsm_ids): Converts a batch of GSM IDs to their corresponding SRX IDs.
    - srx_to_srr_batch(srx_ids): Converts a batch of SRX IDs to their corresponding SRR IDs.
    - main(): Orchestrates the execution of the script.

Error Handling:
    The script includes error handling for HTTP requests. If an error occurs while
    fetching data from the GEO or SRA databases, an exception will be raised with an
    appropriate error message.

Notes:
    - Ensure that the internet connection is stable, as the script fetches data
      from online databases.
    - If no SRR IDs are found for a given GSM ID, the corresponding cell in the
      output will be empty.
    - The script processes each GSM ID sequentially, which may take time depending
      on the number of IDs and network conditions.

Example:
    Input file (samples.csv):
    ```
    GSM
    GSM123456
    GSM123457
    ```

    After running the script, the updated samples.csv might look like this:
    ```
    GSM    SRR
    GSM123456  SRR123456
    GSM123457  SRR123457
    ```

    The output file (SRRs.txt) will contain:
    ```
    ('SRR123456' 'SRR123457')
    """
import argparse
from bs4 import BeautifulSoup
import requests
import pandas as pd
import csv


def gsm_to_srx(gsm_id):
    """
    Fetch the SRX ID associated with a given GSM ID.

    Parameters:
        gsm_id (str): The GSM ID to query.

    Returns:
        str or None: The associated SRX ID or None if not found.
    """
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id}"
    response = requests.get(url)
    response.raise_for_status()  # Raises an error for bad responses
    soup = BeautifulSoup(response.content, "html.parser")
    srx = soup.find("a", string="SRA Run Selector")
    if srx is not None:
        return srx["href"].split("=")[-1]
    return None

def srx_to_srr(srx_id):
    """
    Fetch the SRR ID associated with a given SRX ID.

    Parameters:
        srx_id (str): The SRX ID to query.

    Returns:
        str or None: The associated SRR ID or None if not found.
    """
    url = f"https://www.ncbi.nlm.nih.gov/sra/?term={srx_id}"
    response = requests.get(url)
    response.raise_for_status()  # Raises an error for bad responses
    soup = BeautifulSoup(response.content, "html.parser")
    srr_links = soup.find_all("a", href=True, string=lambda x: x and x.startswith("SRR"))
    if srr_links:
        return srr_links[0].text
    return None

def gsm_to_srx_batch(gsm_ids):
    """
    Convert a batch of GSM IDs to their corresponding SRX IDs.

    Parameters:
        gsm_ids (list): A list of GSM IDs to convert.

    Returns:
        list: A list of SRX IDs corresponding to the input GSM IDs.
    """
    srxs = []
    for gsm in gsm_ids:
        srx = gsm_to_srx(gsm)
        if srx is not None:
            srxs.append(srx)
            print(f"{gsm} -> {srx}")
    return srxs

def srx_to_srr_batch(srx_ids):
    """
    Convert a batch of SRX IDs to their corresponding SRR IDs.

    Parameters:
        srx_ids (list): A list of SRX IDs to convert.

    Returns:
        list: A list of SRR IDs corresponding to the input SRX IDs.
    """
    srrs = []
    for srx in srx_ids:
        srr = srx_to_srr(srx)
        if srr is not None:
            srrs.append(srr)
            print(f"{srx} -> {srr}")
    return srrs


def main():
    # Argument parser for command line arguments
    parser = argparse.ArgumentParser(description="Convert GSM IDs to SRR IDs")
    parser.add_argument("workdir", help='Directory containing sample file with "GSM" column')

    args = parser.parse_args()

    # Paths
    sample_file = args.workdir + "/samples.csv"
    srr_file = args.workdir + "/SRRs.txt"

    # Read GSM IDs from sample file
    df = pd.read_csv(sample_file, sep="\t")
    print(df)

    gsm_ids = df["GSM"].to_list()
    print(f"\nConverting: {len(gsm_ids)} samples")

    # Convert GSM IDs to SRX IDs
    SRXs = gsm_to_srx_batch(gsm_ids)

    # Convert SRX IDs to SRR IDs
    SRRs = srx_to_srr_batch(SRXs)
    print(SRRs)

    df["SRR"] = pd.Series(SRRs)
    df.to_csv(sample_file, sep= "\t", index=False)

    # Write bash-style list of SRR IDs to file
    srr_list = df["SRR"].to_list()
    bash_list = "('" + "' '".join(srr_list) + "')"

    with open(srr_file, 'w', newline='\n') as file:
        writer = csv.writer(file)
        writer.writerow([bash_list])

if __name__ == "__main__":
    main()
