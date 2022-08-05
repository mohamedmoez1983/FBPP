#this module for creating function to make BLAST for predict primer 

# ======================
# imports
# ======================
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO
from textwrap import dedent
import sqlite3






# ======================
# Main function
# ======================


def Primer_BLAST(primer_seq, organism_name):
    # ======================
    # DateBase link 
    # ======================
    # ==================================================== 
    #*******create Database or connect if found******** 
    # ==================================================== 
    conn = sqlite3.connect('vfdb.db' )
    
    
    # ==================================================== 
    #*******Dealing with the database******** 
    # ==================================================== 
    cur = conn.cursor()
    
    # ==================================================== 
    #*******create the table in date base******** 
    # ==================================================== 
    
    
    # ==================================================== 
    # Primer Table
    # ==================================================== 
    cur.execute("""CREATE TABLE  IF NOT EXISTS Primer (
                       primer_seq TEXT PRIMARY KEY,
                       file_name TEXT
                       )""")
    
    
    # ==================================================== 
    #*******Select the data from the database******** 
    # ==================================================== 
    
    cur.execute("SELECT primer_seq FROM  Primer")
    primer_seq_quary = cur.fetchall()
    primer_seq_list = []
    for x in primer_seq_quary:
        primer_seq_list.append(x[0])
        
    cur.execute("""SELECT file_name FROM  Primer WHERE  primer_seq = ?""", (primer_seq,))
    file_name_quary = cur.fetchall()
    if primer_seq in primer_seq_list:
        for blast_qresult in SearchIO.read(file_name_quary[0][0],'blast-xml'):
            if organism_name not in blast_qresult.description:
                return dedent(f"""this primer is not unique
                   the organism match with {blast_qresult.description}
                   the id of blast {blast_qresult.id}""")    
                break
            else:
                return (" this primer is  unique","existed file")
    else:
        result_handle = NCBIWWW.qblast ("blastn", "nt", primer_seq)
        name = str(len(primer_seq_list)+1)+".xml"
        # ==================================================== 
        #*******create Database or connect if found******** 
        # ==================================================== 
        conn = sqlite3.connect('vfdb.db' )
        
        # ==================================================== 
        #*******Dealing with the database******** 
        # ==================================================== 
        cur = conn.cursor()        
        # ====================================================
        # *******Insert the data into the database********
        # ====================================================
        cur.execute('INSERT INTO Primer VALUES (?,?)', (primer_seq, name))   
        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()
        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()        
        with open(name, "w") as out_handle:
            out_handle.write(result_handle.read())
            result_handle.close()          
        for blast_qresult in SearchIO.read(name,'blast-xml'):
            if organism_name not in blast_qresult.description:
                    return dedent(f"""this primer is not unique
                   the organism match with {blast_qresult.description}
                   the id of blast {blast_qresult.id}""") 
                    break
        return (" this primer is  unique","new file")
          
    
# the below lines were created for testing the function                   
#print(Primer_BLAST("TCCTTTGACGGTGCGATGAA", "Salmonella enterica"))     

