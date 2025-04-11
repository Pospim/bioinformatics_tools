from Bio import Entrez
from bs4 import BeautifulSoup
import time, requests, pickle
import pandas as pd
import xml.etree.ElementTree as ET

# Set your email for Entrez
Entrez.email = "marek.pospisil@uochb.cas.cz"

def extract_metadata(soup, tag):
    """
    Extract metadata from the BeautifulSoup object based on the specified tag.

    Parameters:
        soup (BeautifulSoup): The parsed HTML document.
        tag (str): The metadata tag to search for.

    Returns:
        str or None: The extracted metadata value or None if not found.
    """
    field = soup.find("td", string=tag)
    if field is not None:
        return field.find_next_sibling("td").text.strip()
    return None

def fetch_geoids(query):
    """
    Search the GEO database for GEO IDs based on a query.

    Parameters:
        query (str): The search term for GEO database.

    Returns:
        list: A list of GEO IDs found.
    """
    handle = Entrez.esearch(db="gds", term=query, retmax=500)
    record = Entrez.read(handle)
    handle.close()

    geo_ids = record["IdList"]
    print(f"Found {len(geo_ids)} GEO IDs\n", geo_ids)
    return geo_ids

def geo_to_gse(geo_ids):
    """
    Convert GEO IDs to GSE IDs.

    Parameters:
        geo_ids (list): A list of GEO IDs.

    Returns:
        list: A list of corresponding GSE IDs.
    """
    gse_ids = []
    id_string = ",".join(geo_ids)
    handle = Entrez.esummary(db="gds", id=id_string, retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    for record in records:
        gse_ids.append(record["Accession"])
    return gse_ids

def fetch_geo_data(gse_id):
    """
    Fetch detailed GEO data from the web for a given GSE ID.

    Parameters:
        gse_id (str): The GSE ID to fetch data for.

    Returns:
        dict or None: A dictionary containing the GSE data or None if an error occurs.
    """
    try:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.content, "html.parser")
        title = extract_metadata(soup, "Title")
        overall_design = extract_metadata(soup, "Overall design")
        summary = extract_metadata(soup, "Summary")

        print(url)
        if summary is None:
            return None

        return {
            "id": gse_id,
            "title": title,
            "summary": summary,
            "overall_design": overall_design,
            "link": url
        }
    except Exception as e:
        print(f"Error processing GSE ID: {gse_id}")
        print(e)
        return None

def fetch_geo_list(gse_ids):
    """
    Fetch a list containing GEO data (as dict) for a list of GSE IDs.

    Parameters:
        gse_ids (list): A list of GSE IDs.
        include_title (bool): Whether to include the title in the DataFrame.
        include_od (bool): Whether to include the overall design in the DataFrame.
        include_summary (bool): Whether to include the summary in the DataFrame.

    Returns:
        List: List containing the GEO data.
    """
    records = []
    for gse_id in gse_ids:
        try:
            record = fetch_geo_data(gse_id)  # Assume this function returns a dictionary
            if record is not None:
                records.append(record)
                print(f"Processed GSE ID: {gse_id}")
            else:
                print(f"No data for GSE ID: {gse_id}")
        except Exception as e:
            print(f"Error processing GSE ID: {gse_id} - {e}")
        time.sleep(0.3)  # Delay to comply with API guidelines

    if not records:
        raise ValueError("No valid records were fetched.")
    return records

def geo_to_df(records, include_title=True, include_od=True, include_summary=True):
    """
    Convert a list of GEO records to a DataFrame.

    Parameters:
        record_dict (dict): A dictionary containing GEO records.
        include_title (bool): Whether to include the title in the DataFrame.
        include_od (bool): Whether to include the overall design in the DataFrame.
        include_summary (bool): Whether to include the summary in the DataFrame.

    Returns:
        DataFrame: A pandas DataFrame containing the GEO data.
    """
    df = pd.DataFrame(records)

    # Drop columns based on user preferences
    columns_to_drop = []
    if not include_title:
        columns_to_drop.append("title")
    if not include_od:
        columns_to_drop.append("overall_design")
    if not include_summary:
        columns_to_drop.append("summary")

    df.drop(columns=[col for col in columns_to_drop if col in df.columns], inplace=True)


    if "id" in df.columns:
        df.set_index("id", inplace=True)
    else:
        raise ValueError("The 'id' column is missing in the DataFrame.")

    return df

# Example usage - replace with actual data paths
csv_file = "/home/pospim/Desktop/work/bioinformatics/datasets/other/cell_death.csv"
"""
query = '''("SARS-CoV-2"[Title] OR "COVID-19"[Title] OR "coronavirus"[Title])
            AND ("human cells"[All Fields] OR "Homo sapiens"[All Fields])
            AND "Homo sapiens"[Organism] AND ("infection"[All Fields] OR "infected"[All Fields])
            AND ("mRNA"[All Fields]))'''

"""
query = '''("human cells"[All Fields] OR "Homo sapiens"[All Fields])
            AND "Homo sapiens"[Organism] AND "mRNA"[All Fields]
            AND ("apoptosis"[All Fields] OR "necrosis"[All Fields] OR
            "cell death"[All Fields])
            '''

geo_ids = fetch_geoids(query)
gse_ids = geo_to_gse(geo_ids)
print(len(geo_ids), len(gse_ids))

records = fetch_geo_list(gse_ids=gse_ids)


df = geo_to_df(records)
print(df)

df.to_csv(csv_file, sep=",", encoding="utf-8")
