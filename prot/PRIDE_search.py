import requests
import os
import pandas as pd
import glob
import argparse

def get_project_info(data):
    project_descriptions = []

    if isinstance(data,list):
        projects = data
    else:
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

def main():
    parser = argparse.ArgumentParser(description="Search PRIDE database for proteomics projects")
    parser.add_argument("--query", required=True, help="Search query string (e.g., 'influenza')")
    parser.add_argument("--output_dir", required=True, help="Output directory for resulting dataframe.")
    parser.add_argument("--max_pages", type=int, default=5, help="Page number for results (default: 0)")
    parser.add_argument("--page_size", type=int, default=100, help="Number of results per page (default: 100)")
    parser.add_argument("--no-title", action="store_true", help="Exclude title column in output")
    parser.add_argument("--no-desc", action="store_true", help="Exclude description column in output")
    parser.add_argument("--no-spp", action="store_true", help="Exclude sample processing protocol column in output")
    parser.add_argument("--no-dpp", action="store_true", help="Exclude data processing protocol column in output")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    out_file = os.path.join(args.output_dir, f"{args.query}.csv")
    all_results = []
    print(f"Searching PRIDE database for '{args.query}'...")

    for page in range(0, args.max_pages):
        print(f"Fetching page {page}...")
        projects = fetch_pride_df(
            keyword=args.query,
            include_title=not args.no_title,
            include_desc=not args.no_desc,
            include_spp=not args.no_spp,
            include_dpp=not args.no_dpp,
            page=page,
            page_size=args.page_size
        )
        if projects is not None and not projects.empty:
            print(f"Found {len(projects)} projects on page {page}")
            all_results.append(projects)
        else:
            print(f"No more results found on page {page}")
            break
    if all_results:
        combined_df = pd.concat(all_results)
        combined_df = combined_df.loc[~combined_df.index.duplicated(keep='first')]
        print(f"Total projects found: {len(combined_df)}")

        combined_df.to_csv(out_file, sep=",", encoding="utf-8")
        print(f"Resulting dataframe saved => {out_file}")
    else:
        print(f"No results found for '{args.query}'")

if __name__ == "__main__":
    main()

"""
# Example usage - replace with actual data paths
keyword = "influenza"
page = 0
csv_file = f"~/Desktop/work/bioinformatics/tools/proteomics/{keyword}_{page}.csv"

projects = fetch_pride_df(keyword, include_title=True, include_desc=True, include_spp=True, include_dpp=True, page=page, page_size=100)
print(projects)

projects.to_csv(csv_file, sep=",", encoding="utf-8")
"""
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