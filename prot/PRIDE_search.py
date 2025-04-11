import requests
import os
import pandas as pd
import glob

def get_project_info(data):
    project_descriptions = []

    projects = data.get('_embedded', {}).get('compactprojects', [])
    if not projects:
        print("No projects found.")
        return None
    for project in projects:
        project_data = {
            "accession": project.get("accession", "N/A"),
            "title": project.get("title", "N/A"),
            "description": project.get("projectDescription", "N/A"),
            "sample_processing_protocol": project.get("sampleProcessingProtocol", "N/A"),
            "data_processing_protocol": project.get("dataProcessingProtocol", "N/A"),
        }
        project_descriptions.append(project_data)
    return project_descriptions

def fetch_pride_projects(keyword, page, page_size):
    """
    Fetches project descriptions from PRIDE database based on a keyword.

    Parameters:
    - keyword (str): The keyword for the PRIDE database search.
    - page (int): The page number of results to fetch.
    - page_size (int): The number of projects per page.

    Returns:
    - List of dictionaries with project titles and descriptions.
    """
    url = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects"

    headers = {"Accept": "application/json"}

    params = {
        "keyword": keyword,
        "filter": "organisms=Homo sapiens (human)",
        "page": page,
        "pageSize": page_size,
        "sortDirection": "DESC"
    }
    response = requests.get(url, params=params, headers=headers)

    if response.status_code == 200:
        return response.json()
    else:
        print("Error fetching data from PRIDE.")
        return None

def fetch_pride_df(keyword, include_title=True, include_desc=True, include_spp=True, include_dpp=True, page=0, page_size=100):

    records = fetch_pride_projects(keyword, page, page_size)
    if not records:
        return None

    descriptions = get_project_info(records)

    df = pd.DataFrame(descriptions)
    # Drop columns based on user preferences
    if not include_title:
        df.drop(columns=["title"], inplace=True)
    if not include_desc:
        df.drop(columns=["description"], inplace=True)
    if not include_spp:
        df.drop(columns=["sample_processing_protocol"], inplace=True)
    if not include_dpp:
        df.drop(columns=["data_processing_protocol"], inplace=True)

    df.set_index("accession", inplace=True)
    return df

# Example usage - replace with actual data paths
keyword = "influenza,flu,early"
page = 0
csv_file = f"~/Desktop/work/bioinformatics/proteomics/{keyword}_{page}.csv"

projects = fetch_pride_df(keyword, include_title=True, include_desc=True, include_spp=True, include_dpp=True, page=page, page_size=100)
print(projects)

projects.to_csv(csv_file, sep="\t", encoding="utf-8")
"""
dir = os.path.expanduser("~/Desktop/work/bioinformatics/proteomics/")
file_pattern = "sars-cov-2_*.csv"
csv_files = glob.glob(os.path.join(dir, file_pattern))

dfs = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file, sep="\t")
    dfs.append(df)
print(dfs)
combined_df = pd.concat(dfs).drop_duplicates()
combined_df.to_csv(f"{dir}sars-cov-2_combined.csv", sep="\t", encoding="utf-8")
"""