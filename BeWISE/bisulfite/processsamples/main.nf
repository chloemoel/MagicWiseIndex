#!/usr/bin/env nextflow

process PROCESS_METHYLATION_BIS {
    publishDir params.outdir, mode: 'copy'
    
    input:
        val batch_correction
        path bisulfite_files
        path additional_data
        
        
    output:
        path "m_values_processed.csv", emit: csv

    script:
        def batch = batch_correction.collect { batch -> "\"${batch}\"" }.join(", ")
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np
        from inmoose.pycombat import pycombat_norm
        import sqlite3

        conn = sqlite3.connect("bisulfite.db")
        cur = conn.cursor()

        cur.execute('''
        CREATE TABLE IF NOT EXISTS methylation (
            chr TEXT,
            start INTEGER,
            stop INTEGER,
            PRIMARY KEY (chr, start, stop)
        )
        ''')
        conn.commit()
        
        def add_sample(file, conn):
            df = pd.read_csv(file, nrows=1)
            sample_name = df.columns[-1]

            cur = conn.cursor()

            # Add column if not exists
            cur.execute(f"ALTER TABLE methylation ADD COLUMN '{sample_name}' REAL")
            conn.commit()

            for chunk in pd.read_csv(file, header=None, skiprows=1, chunksize=100000):
                chunk.columns = ["chr", "start", "stop", sample_name]

                # Insert rows if they don't exist
                rows = chunk[["chr", "start", "stop"]].itertuples(index=False, name=None)
                cur.executemany('''
                    INSERT OR IGNORE INTO methylation (chr, start, stop)
                    VALUES (?, ?, ?)
                ''', rows)

                # Update values
                rows = chunk.itertuples(index=False, name=None)
                cur.executemany(f'''
                    UPDATE methylation
                    SET '{sample_name}' = ?
                    WHERE chr = ? AND start = ? AND stop = ?
                ''', [(r[3], r[0], r[1], r[2]) for r in rows])

                conn.commit()

        files = "${bisulfite_files}".split()
        for file in files:
            add_sample(file, conn)


        m_values = pd.read_sql_query("SELECT * FROM methylation", conn)
        conn.close()

        # Filter out sites with more than 10% missing data and fill NaN with mean of each probe
        thresh = int(len(m_values.iloc[:, 3:].columns) * 0.1)
        m_values.dropna(thresh=thresh, axis=1, inplace=True)
        m_values.iloc[:, 3:] = m_values.iloc[:, 3:].apply(lambda row: row.fillna(row.mean()),axis=1)



        if "${batch_correction}" != "[]":
            if "${additional_data}" != "":
                additional_data = pd.read_csv("${additional_data}", index_col=0, dtype={0:str})
        
            # Merge sample info and m_values so m values are in right order 
            m_and_info = m_values.T.merge(additional_data, left_index=True, right_index=True)
            header = m_and_info.index.astype(object)
            m_values = m_and_info.drop(columns=additional_data.columns).T 

            # Perform batch correction for each batch in the list
            for b in [${batch}]:
                m_values = pycombat_norm(m_values, m_and_info[b], na_cov_action="remove")           

       
        if len(m_values.columns) != len(set(m_values.columns)):
            print("WARNING! multiple samples found for one or more subjects. We will take the average of their score")
            m_values_ = m_values.T
            m_values_["Header"] = header
            m_values = m_values_.groupby(by="Header",group_keys=False, as_index=True, sort=False).mean().T

        # Save the processed m_values with samples as columns and probes as rows
        m_values.to_csv("m_values_processed.csv", index=False)
        """
}
