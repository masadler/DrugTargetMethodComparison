from urllib.request import Request, urlopen
import pandas as pd
from bs4 import BeautifulSoup
import lxml

df = pd.read_csv("data/GWAS_indication.csv")
df = (df.set_index(df.columns.drop('DBCOND_ID',1).tolist())
                .DBCOND_ID.str.split('; ', expand=True)
                .stack()
                .reset_index()
                .rename(columns={0:'DBCOND_ID'})
                .loc[:, df.columns]
                )

diseases = list(df.Acronym.unique())

for d in diseases:

    print(d)

    indications = list(df["DBCOND_ID"][df.Acronym == d])

    db_id = []
    db_name = []

    for ind in indications:

        url = 'https://www.drugbank.ca/indications/'+ ind
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        page = urlopen(req).read()
        #page = urllib.request.urlopen(url)

        soup = BeautifulSoup(page, features="lxml")

        table = soup.find('table', class_= "table table-condensed table-striped datatable").find('tbody')

        for row in table.find_all('tr'):
            cells = row.find_all('td')
            if len(cells) == 3:
                db_id.append(cells[0].string)
                db_name.append(cells[1].string)

    result_df = pd.DataFrame({"drugbank_id": db_id, "db_name": db_name})
    result_df = result_df.drop_duplicates()

    print(result_df.head())

    result_df.to_csv("output/drug_targets/" + d + "_drugs_drugbank.csv", index = False)