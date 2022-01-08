# ======================
# imports
# ======================

from tkinter import *
from tkinter import ttk
from tkinter import scrolledtext
from tkinter.scrolledtext import ScrolledText
import PIL
from PIL import Image, ImageTk
import pydna
from pydna.dseqrecord import Dseqrecord
from pydna.readers import read
from pydna.amplify import pcr
from Bio.Seq import Seq
import primer3
from pydna.gel import weight_standard_sample, Gel1
#from matplotlib import savefig
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from numpy import arange, sin, pi
from Primer_Blast import Primer_BLAST
import sqlite3
from Bio import SeqIO
from tkinter import filedialog
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
st = weight_standard_sample('1kb+_GeneRuler')
#the below line for testing in the end project will connnected to datebase


class Parent(Tk):
    # ======================
    # Create instance
    # ======================
    def __init__(self):
        super(Parent,self).__init__()
        self.geometry('675x690+250+25')
        # =========================
        # Disable resizing the GUI
        # =========================
        self.resizable(0, 0)
        # ======================
        # Add a title
        # ======================
        self.title('Food Primer / Probe Predict System')
        # ======================
        # Add a Tab
        # ======================
        self.tabControl = ttk.Notebook(self)
        self.tab1 = ttk.Frame(self.tabControl)
        self.tabControl.add(self.tab1, text="Add New Gene")
        self.tab2 = ttk.Frame(self.tabControl)
        self.tabControl.add(self.tab2, text="Select primer")
        self.tabControl.grid(column=0, row=0)
        self.Add_image_tab1()
        self.Selection_quary()
        self.Entry_Parmeter()


    # =========================
    # Create Database function
    # =========================
    def open_file(self):
        self.label_Accession_name.configure(text="")
        self.label_protein_id_name.configure(text="")
        self.label_protein_name_name.configure(text="")
        self.label_Locus_name.configure(text="")
        self.label_DBlink_name.configure(text="")
        self.label_DBlink_name.configure(text="")
        self.label_Taxon_name.configure(text="")
        self.label_Version_name.configure(text="")
        self.label_Gene_id_name.configure(text="")
        self.label_Gene_name_name.configure(text="")
        self.label_Gene_region_name.configure(text="")
        self.label_Gene_Size_name.configure(text="")
        self.label_Genus_name.configure(text="")
        self.label_Serovar_name.configure(text="")
        self.label_Species_name.configure(text="")
        self.label_Sub_Species_name.configure(text="")
        self.filename = filedialog.askopenfile(initialdir="/C", title="Select A Genbank File", filetypes=(("Genbank files", "*.gb"), ("All file", "*.*")))
        record = SeqIO.read(self.filename, "genbank")

        # ====================================================
        # date fill into the database
        # ====================================================

        # ====================================================
        # 1- for finding "locus" or "Accession " in All table
        # ====================================================
        # print(record.id)
        Accession_no = record.annotations['accessions'][0]

        # ====================================================
        # 2- for finding "Genus name " in organis_genus table
        # ====================================================
        # print(record.annotations['source'])
        Genus_name = str(record.annotations['organism'].split()[0])

        # ====================================================
        # 3 - for finding "Species name " in organis_Species table
        # ====================================================
        # print(record.annotations['source'])
        Species_name = str(record.annotations['organism'].split()[1])

        # ====================================================
        # 4- for finding "sub species" in organis table
        # ====================================================
        for feature in record.features:
            if 'sub_species' in feature.qualifiers:
                Sup_Species = feature.qualifiers.get("sub_species")[0]
                break
            else:
                Sup_Species = "None"

        # ====================================================
        # 5- for finding "serovar" in organis table
        # ====================================================
        for feature in record.features:
            if 'serovar' in feature.qualifiers:
                Serovar = feature.qualifiers.get("serovar")[0]
                break
            else:
                Serovar = "None"

        # ====================================================
        # 6- for finding "strain" in organis table
        # ====================================================
        for feature in record.features:
            if 'strain' in feature.qualifiers:
                Strain = feature.qualifiers.get("strain")[0]
                break
            else:
                Strain = "None"

        # ====================================================
        # 7- for finding "molecule type" in organis table
        # ====================================================
        for feature in record.features:
            if 'mol_type' in feature.qualifiers:
                Molecule_type = feature.qualifiers.get("mol_type")[0]
                break

            else:
                Molecule_type = "None"

            # ====================================================
        # 8- for finding "DBLINK" in organis table
        # ====================================================
        DBLINK = str(record.dbxrefs)

        # ====================================================
        # 9- for finding "Organism Database Ref" in organis table
        # ====================================================
        for feature in record.features:
            if 'strain' in feature.qualifiers:
                taxon = feature.qualifiers.get("db_xref")[0]
                break
            else:
                taxon = "None"

        # ====================================================
        # 10- for finding "version " in organis table
        # ====================================================
        Version = record.annotations['sequence_version']

        # ====================================================
        # 1- for finding "locus_tag" in Gene table
        # ====================================================
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers.get("locus_tag")[0]
                break
            else:
                locus_tag = "None"

        # ====================================================
        # 3- for finding "Gene ID" in Gene table
        # ====================================================
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                if 'db_xref' in feature.qualifiers:
                    Gene_ID = feature.qualifiers.get("db_xref")[0]
                    break
                else:
                    Gene_ID = "None"

        # ====================================================
        # 4- for finding "Size" in Gene table
        # ====================================================
        Size = len(record.seq)

        # ====================================================
        # 5- for finding "Region" in Gene table
        # ====================================================
        Region = record.annotations['accessions'][2]

        # ====================================================
        # 6 & 7- for finding "Product Name & Protein ID" in Gene table
        # ====================================================

        for feature in record.features:
            if 'product' in feature.qualifiers:
                Prodcut_Name = feature.qualifiers.get("product")[0]

        for feature in record.features:
            if 'protein_id' in feature.qualifiers:
                Protein_ID = feature.qualifiers.get("protein_id")[0]
                break
            else:
                Protein_ID = "None"

        # ====================================================
        # 2- for finding "Gene name" in Gene table
        # ====================================================
        for feature in record.features:
            if 'gene' in feature.qualifiers:
                Gene_Name = feature.qualifiers.get("gene")[0]
                break
            else:
                Gene_Name = Prodcut_Name + " gene"

        # ====================================================
        # 9- for finding "Sequence" in Gene table
        # ====================================================
        Sequence = str(record.seq)

        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')

        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()

        # ====================================================
        # *******create the table in date base********
        # ====================================================

        # ====================================================
        # Organism Table
        # ====================================================
        self.cur.execute("""CREATE TABLE  IF NOT EXISTS Organism (
                           Accession TEXT PRIMARY KEY,
                           Genus TEXT,
                           Species TEXT,
                           Sub_Species TEXT,
                           Serovar TEXT,
                           Strain TEXT,
                           Molecule_type TEXT,
                           DBlink TEXT,
                           Taxon TEXT , 
                           Version INTEGER
                           )""")

        # ====================================================
        # Gene Table
        # ====================================================
        self.cur.execute("""CREATE TABLE IF NOT EXISTS Gene (
                           locus_tag TEXT PRIMARY KEY,
                           Gene_Name TEXT, 
                           Gene_ID TEXT,
                           Size INTEGER,
                           Region TEXT,
                           Prodcut_Name TEXT,
                           Protein_ID TEXT,
                           Accession TEXT,
                           Sequence TEXT
                           )""")

        # ====================================================
        # *******Insert the data into the database********
        # ====================================================
        # if Accession_no  in cur.execute('''SELECT Accession FROM Organism '''):
        try:
            self.cur.execute('INSERT INTO Organism VALUES (?,?,?,?,?,?,?,?,?,?)', (
                Accession_no, Genus_name, Species_name, Sup_Species, Serovar, Strain, Molecule_type, DBLINK, taxon,
                Version))
        except:
            pass
        try:
            self.cur.execute('INSERT INTO Gene VALUES (?,?,?,?,?,?,?,?,?)',
                        (locus_tag, Gene_Name, Gene_ID, Size, Region, Prodcut_Name, Protein_ID, Accession_no,
                         Sequence))
        except:
            pass


        # ========================================================
        # label_configure for Database in  Add gene Frame
        # ========================================================
        self.label_Accession_name.configure(text=Accession_no)
        self.label_Version_name.configure(text=Version)
        self.label_Genus_name.configure(text=Genus_name)
        self.label_Species_name.configure(text=Species_name)
        self.label_Sub_Species_name.configure(text=Sup_Species)
        self.label_Serovar_name.configure(text=Serovar)
        self.label_Taxon_name.configure(text=taxon)
        self.label_DBlink_name.configure(text=DBLINK)
        self.label_Locus_name.configure(text=locus_tag)
        self.label_Gene_id_name.configure(text=Gene_ID)
        self.label_Gene_name_name.configure(text=Gene_Name)
        self.label_Gene_Size_name.configure(text=Size)
        self.label_Gene_region_name.configure(text=Region)
        self.label_protein_name_name.configure(text=Prodcut_Name)
        self.label_protein_id_name.configure(text=Protein_ID)

        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()

        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()

    # =========================
    # Create Database function
    # =========================
    def retrive_gene(self):
        self.label_Accession_name.configure(text="")
        self.label_protein_id_name.configure(text="")
        self.label_protein_name_name.configure(text="")
        self.label_Locus_name.configure(text="")
        self.label_DBlink_name.configure(text="")
        self.label_DBlink_name.configure(text="")
        self.label_Taxon_name.configure(text="")
        self.label_Version_name.configure(text="")
        self.label_Gene_id_name.configure(text="")
        self.label_Gene_name_name.configure(text="")
        self.label_Gene_region_name.configure(text="")
        self.label_Gene_Size_name.configure(text="")
        self.label_Genus_name.configure(text="")
        self.label_Serovar_name.configure(text="")
        self.label_Species_name.configure(text="")
        self.label_Sub_Species_name.configure(text="")
        Ent_Accession = self.entery_Accession.get()
        Entrez.email = "m.moez.1983@gmail.com"  # Always tell NCBI who you are
        handle = Entrez.efetch(db="nucleotide", id= Ent_Accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")

        # ====================================================
        # date fill into the database
        # ====================================================

        # ====================================================
        # 1- for finding "locus" or "Accession " in All table
        # ====================================================
        # print(record.id)
        Accession_no = record.annotations['accessions'][0]

        # ====================================================
        # 2- for finding "Genus name " in organis_genus table
        # ====================================================
        # print(record.annotations['source'])
        Genus_name = str(record.annotations['organism'].split()[0])

        # ====================================================
        # 3 - for finding "Species name " in organis_Species table
        # ====================================================
        # print(record.annotations['source'])
        Species_name = str(record.annotations['organism'].split()[1])

        # ====================================================
        # 4- for finding "sub species" in organis table
        # ====================================================
        for feature in record.features:
            if 'sub_species' in feature.qualifiers:
                Sup_Species = feature.qualifiers.get("sub_species")[0]
                break
            else:
                Sup_Species = "None"

        # ====================================================
        # 5- for finding "serovar" in organis table
        # ====================================================
        for feature in record.features:
            if 'serovar' in feature.qualifiers:
                Serovar = feature.qualifiers.get("serovar")[0]
                break
            else:
                Serovar = "None"

        # ====================================================
        # 6- for finding "strain" in organis table
        # ====================================================
        for feature in record.features:
            if 'strain' in feature.qualifiers:
                Strain = feature.qualifiers.get("strain")[0]
                break
            else:
                Strain = "None"

        # ====================================================
        # 7- for finding "molecule type" in organis table
        # ====================================================
        for feature in record.features:
            if 'mol_type' in feature.qualifiers:
                Molecule_type = feature.qualifiers.get("mol_type")[0]
                break

            else:
                Molecule_type = "None"

            # ====================================================
        # 8- for finding "DBLINK" in organis table
        # ====================================================
        DBLINK = str(record.dbxrefs)

        # ====================================================
        # 9- for finding "Organism Database Ref" in organis table
        # ====================================================
        for feature in record.features:
            if 'strain' in feature.qualifiers:
                taxon = feature.qualifiers.get("db_xref")[0]
                break
            else:
                taxon = "None"

        # ====================================================
        # 10- for finding "version " in organis table
        # ====================================================
        Version = record.annotations['sequence_version']

        # ====================================================
        # 1- for finding "locus_tag" in Gene table
        # ====================================================
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers.get("locus_tag")[0]
                break
            else:
                locus_tag = "None"

        # ====================================================
        # 3- for finding "Gene ID" in Gene table
        # ====================================================
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                if 'db_xref' in feature.qualifiers:
                    Gene_ID = feature.qualifiers.get("db_xref")[0]
                    break
                else:
                    Gene_ID = "None"
            elif locus_tag  == "None":
                Gene_ID = "None"

        # ====================================================
        # 4- for finding "Size" in Gene table
        # ====================================================
        Size = len(record.seq)

        # ====================================================
        # 5- for finding "Region" in Gene table
        # ====================================================
        for annotaed in record.annotations:
            if "REGION:" in record.annotations:
                Region = record.annotations['accessions'][1]
                break
            else:
                Region = "None"
                break

        # ====================================================
        # 6 & 7- for finding "Product Name & Protein ID" in Gene table
        # ====================================================

        for feature in record.features:
            if 'product' in feature.qualifiers:
                Prodcut_Name = feature.qualifiers.get("product")[0]

        for feature in record.features:
            if 'protein_id' in feature.qualifiers:
                Protein_ID = feature.qualifiers.get("protein_id")[0]
                break
            else:
                Protein_ID = "None"

        # ====================================================
        # 2- for finding "Gene name" in Gene table
        # ====================================================
        for feature in record.features:
            if 'gene' in feature.qualifiers:
                Gene_Name = feature.qualifiers.get("gene")[0]
                break
            else:
                Gene_Name = Prodcut_Name + " gene"

        # ====================================================
        # 9- for finding "Sequence" in Gene table
        # ====================================================
        Sequence = str(record.seq)

        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')

        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()

        # ====================================================
        # *******create the table in date base********
        # ====================================================

        # ====================================================
        # Organism Table
        # ====================================================
        self.cur.execute("""CREATE TABLE  IF NOT EXISTS Organism (
                           Accession TEXT PRIMARY KEY,
                           Genus TEXT,
                           Species TEXT,
                           Sub_Species TEXT,
                           Serovar TEXT,
                           Strain TEXT,
                           Molecule_type TEXT,
                           DBlink TEXT,
                           Taxon TEXT , 
                           Version INTEGER
                           )""")

        # ====================================================
        # Gene Table
        # ====================================================
        self.cur.execute("""CREATE TABLE IF NOT EXISTS Gene (
                           locus_tag TEXT PRIMARY KEY,
                           Gene_Name TEXT, 
                           Gene_ID TEXT,
                           Size INTEGER,
                           Region TEXT,
                           Prodcut_Name TEXT,
                           Protein_ID TEXT,
                           Accession TEXT,
                           Sequence TEXT
                           )""")

        # ====================================================
        # *******Insert the data into the database********
        # ====================================================
        # if Accession_no  in cur.execute('''SELECT Accession FROM Organism '''):
        try:
            self.cur.execute('INSERT INTO Organism VALUES (?,?,?,?,?,?,?,?,?,?)', (
                Accession_no, Genus_name, Species_name, Sup_Species, Serovar, Strain, Molecule_type, DBLINK, taxon,
                Version))
        except:
            pass
        try:
            self.cur.execute('INSERT INTO Gene VALUES (?,?,?,?,?,?,?,?,?)',
                             (locus_tag, Gene_Name, Gene_ID, Size, Region, Prodcut_Name, Protein_ID, Accession_no,
                              Sequence))
        except:
            pass

        # ========================================================
        # label_configure for Database in  Add gene Frame
        # ========================================================
        self.label_Accession_name.configure(text=Accession_no)
        self.label_Version_name.configure(text=Version)
        self.label_Genus_name.configure(text=Genus_name)
        self.label_Species_name.configure(text=Species_name)
        self.label_Sub_Species_name.configure(text=Sup_Species)
        self.label_Serovar_name.configure(text=Serovar)
        self.label_Taxon_name.configure(text=taxon)
        self.label_DBlink_name.configure(text=DBLINK)
        self.label_Locus_name.configure(text=locus_tag)
        #self.label_Gene_id_name.configure(text=Gene_ID)
        self.label_Gene_name_name.configure(text=Gene_Name)
        self.label_Gene_Size_name.configure(text=Size)
        #self.label_Gene_region_name.configure(text=Region)
        self.label_protein_name_name.configure(text=Prodcut_Name)
        self.label_protein_id_name.configure(text=Protein_ID)

        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()

        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()

    # =========================
    # 1- Create  function
    # =========================

    def rad_gene_Call(self):
        rad_gene_Sel = self.radvar_gene.get()
        if rad_gene_Sel == 1:
            self.label_retrieve_button.configure(state='normal')
            self.entery_Accession.configure(state='normal')
            self.label_open_file_button.configure(state='disable')
            self.label_Accession_name.configure(text = "")
            self.label_protein_id_name.configure(text = "")
            self.label_protein_name_name.configure(text = "")
            self.label_Locus_name.configure(text = "")
            self.label_DBlink_name.configure(text = "")
            self.label_DBlink_name.configure(text = "")
            self.label_Taxon_name.configure(text = "")
            self.label_Version_name.configure(text = "")
            self.label_Gene_id_name.configure(text = "")
            self.label_Gene_name_name.configure(text = "")
            self.label_Gene_region_name.configure(text = "")
            self.label_Gene_Size_name.configure(text = "")
            self.label_Genus_name.configure(text = "")
            self.label_Serovar_name.configure(text = "")
            self.label_Species_name.configure(text = "")
            self.label_Sub_Species_name.configure(text = "")

        elif rad_gene_Sel == 2:
            self.label_open_file_button.configure(state='normal')
            self.label_retrieve_button.configure(state='disable')
            self.entery_Accession.delete(0, END)
            self.entery_Accession.configure(state='disable')
            self.label_Accession_name.configure(text="")
            self.label_protein_id_name.configure(text="")
            self.label_protein_name_name.configure(text="")
            self.label_Locus_name.configure(text="")
            self.label_DBlink_name.configure(text="")
            self.label_DBlink_name.configure(text="")
            self.label_Taxon_name.configure(text="")
            self.label_Version_name.configure(text="")
            self.label_Gene_id_name.configure(text="")
            self.label_Gene_name_name.configure(text="")
            self.label_Gene_region_name.configure(text="")
            self.label_Gene_Size_name.configure(text="")
            self.label_Genus_name.configure(text="")
            self.label_Serovar_name.configure(text="")
            self.label_Species_name.configure(text="")
            self.label_Sub_Species_name.configure(text="")
        else:
            self.label_open_file_button.configure(state='disable')
            self.label_retrieve_button.configure(state='disable')



    def radCall(self):
        self.radSel = self.radvar_Sel.get()
        if self.radSel == 1:
            self.c_left_primer.configure(state='normal')
            self.c_right_primer.configure(state='normal')
            self.c_probe.configure(state='normal')
        elif self.radSel == 2 or self.radSel == 3:
            self.c_left_primer.configure(state='disabled')
            self.c_right_primer.configure(state='disabled')
            self.c_probe.configure(state='disabled')

    def probe_Call(self):
        self.rad_probe_Sel = self.var_p.get()
        if  self.rad_probe_Sel == 1:
            self.c_left_primer.configure(state='disabled')
            self.c_left_primer.deselect()
            self.c_right_primer.configure(state='disabled')
            self.c_right_primer.deselect()
            self.entery_Size_1.configure(state='disabled')
            self.entery_Size_2.configure(state='disabled')
            self.entery_Size_3.configure(state='disabled')
            self.entery_GC_1.configure(state='disabled')
            self.entery_GC_2.configure(state='disabled')
            self.entery_Max_Pair_Compl_1.configure(state='disabled')
            self.entery_Max_Pair_Compl_2.configure(state='disabled')
            self.entery_Max_Self_Compl_1.configure(state='disabled')
            self.entery_Max_Self_Compl_2.configure(state='disabled')
            self.entery_Poly_N.configure(state='disabled')
            self.entery_Poly_x.configure(state='disabled')
            self.entery_con_di_cat.configure(state='disabled')
            self.entery_con_mon_cat.configure(state='disabled')
            self.entery_con_dNTP.configure(state='disabled')
            self.entery_con_AO.configure(state='disabled')
            self.entery_Temp_1.configure(state='disabled')
            self.entery_Temp_2.configure(state='disabled')
            self.entery_Temp_3.configure(state='disabled')
            self.label_Product_Size.configure(text = "Probe Size")
            self.label_Primer_return.configure(text="# of Probe to return")
        else:
            self.c_left_primer.configure(state='normal')
            self.c_right_primer.configure(state='normal')
            self.entery_Size_1.configure(state='normal')
            self.entery_Size_2.configure(state='normal')
            self.entery_Size_3.configure(state='normal')
            self.entery_GC_1.configure(state='normal')
            self.entery_GC_2.configure(state='normal')
            self.entery_Max_Pair_Compl_1.configure(state='normal')
            self.entery_Max_Pair_Compl_2.configure(state='normal')
            self.entery_Max_Self_Compl_1.configure(state='normal')
            self.entery_Max_Self_Compl_2.configure(state='normal')
            self.entery_Poly_N.configure(state='normal')
            self.entery_Poly_x.configure(state='normal')
            self.entery_con_di_cat.configure(state='normal')
            self.entery_con_mon_cat.configure(state='normal')
            self.entery_con_dNTP.configure(state='normal')
            self.entery_con_AO.configure(state='normal')
            self.entery_Temp_1.configure(state='normal')
            self.entery_Temp_2.configure(state='normal')
            self.entery_Temp_3.configure(state='normal')
            self.label_Product_Size.configure(text="PCR Product Size")
            self.label_Primer_return.configure(text="# of primers to return")
        # ======================
        # Add image to Tab 1
        # ======================
    def Add_image_tab1(self):
        self.image_1 = Image.open("E:/mohamed/PHD/bioinformatics/power point/Gene1.jpg")
        self.photo_1 = ImageTk.PhotoImage(self.image_1)
        self.lab_1 = Label(self.tab1, image=self.photo_1)
        self.lab_1.grid(column=0, row=2)
        # ==========================================
        # Add Database Frame to Tab 1
        # ==========================================
        self.label_open_file_frame = LabelFrame(self.tab1, text="Fetch or Open A Genbank File", fg="red")
        self.label_open_file_frame.place(x=10, y=150, width=650, height=500)

        # ====================================================================
        #  Add "Fetch  / open  Radiobutton" to open file frame
        # ====================================================================
        self.radvar_gene = IntVar()
        R_fetch_genbank = Radiobutton(self.label_open_file_frame, text=" Fetch GenBank File from NCBI through Accession No.", variable=self.radvar_gene, value=1, command=self.rad_gene_Call)
        R_fetch_genbank.grid(column=0, row=0, sticky="W")
        R_open_file = Radiobutton(self.label_open_file_frame, text=" Upload Existing GenBank File", variable=self.radvar_gene, value=2, command=self.rad_gene_Call)
        R_open_file.grid(column=0, row=1, sticky="W")

        # ====================================================================
        #  Add "Entery of Accession No." to open file frame
        # ====================================================================
        self.e_Acession = StringVar()
        self.entery_Accession = Entry(self.label_open_file_frame, width=25, textvariable=self.e_Acession)
        self.entery_Accession.grid(column=2, row=0)
        self.e_Acession.set("Entry of Accession No.")
        self.entery_Accession.configure(state='disable')


        # =========================================
        # Add a "Genbank file information" Frame to Tab 1
        # =========================================
        self.label_frame_information = LabelFrame(self.label_open_file_frame, text="Genbank file information", fg="red")
        self.label_frame_information.place(x=10, y=55, width=600, height=420)

        # =========================================================================
        # Adding a "Genbank information " Label, and data to"Genbank file information"  Frame
        # =========================================================================
        # =========================================================================
        # Adding a "Organism Information " Frame to"Genbank file information"  Frame
        # =========================================================================
        self.label_organism = LabelFrame(self.label_frame_information, text="Organism Information ", fg="blue")
        self.label_organism.place(x=10, y=10, width=500, height=190)
        # =========================================================================
        # Adding a "Organism Information " Label, and data to"Organism information"  Frame
        # =========================================================================
        self.label_Accession = Label(self.label_organism, text="Accession: ")
        self.label_Accession.grid(column=0, row=1, sticky="W")
        self.label_Accession_name = Label(self.label_organism, text="")
        self.label_Accession_name.grid(column=1, row=1, sticky="W")
        self.label_Version = Label(self.label_organism, text="Version: ")
        self.label_Version.grid(column=0, row=2, sticky="W")
        self.label_Version_name = Label(self.label_organism, text="")
        self.label_Version_name.grid(column=1, row=2, sticky="W")
        self.label_Genus = Label(self.label_organism, text="Genus: ")
        self.label_Genus.grid(column=0, row=3, sticky="W")
        self.label_Genus_name = Label(self.label_organism, text="")
        self.label_Genus_name.grid(column=1, row=3, sticky="W")
        self.label_Species = Label(self.label_organism, text="Species: ")
        self.label_Species.grid(column=0, row=4, sticky="W")
        self.label_Species_name = Label(self.label_organism, text="")
        self.label_Species_name.grid(column=1, row=4, sticky="W")
        self.label_Sub_Species = Label(self.label_organism, text="Sub_Species: ")
        self.label_Sub_Species.grid(column=0, row=5, sticky="W")
        self.label_Sub_Species_name = Label(self.label_organism, text="")
        self.label_Sub_Species_name.grid(column=1, row=5, sticky="W")
        self.label_Serovar = Label(self.label_organism, text="Serovar: ")
        self.label_Serovar.grid(column=0, row=6, sticky="W")
        self.label_Serovar_name = Label(self.label_organism, text="")
        self.label_Serovar_name.grid(column=1, row=6, sticky="W")
        self.label_Taxon = Label(self.label_organism, text="Taxon: ")
        self.label_Taxon.grid(column=0, row=7, sticky="W")
        self.label_Taxon_name = Label(self.label_organism, text="")
        self.label_Taxon_name.grid(column=1, row=7, sticky="W")
        self.label_DBlink = Label(self.label_organism, text="DBlink: ")
        self.label_DBlink.grid(column=0, row=8, sticky="W")
        self.label_DBlink_name = Label(self.label_organism, text="")
        self.label_DBlink_name.grid(column=1, row=8, sticky="W")

        # =========================================================================
        # Adding a "Gene Information " Frame to"Genbank file information"  Frame
        # =========================================================================
        self.label_Gene = LabelFrame(self.label_frame_information, text="Gene Information ", fg="blue")
        self.label_Gene.place(x=10, y=220, width=500, height=170)
        # =========================================================================
        # Adding a "Gene Information " Label, and data to"Gene information"  Frame
        # =========================================================================
        self.label_Locus = Label(self.label_Gene, text="Locus_tag: ")
        self.label_Locus.grid(column=0, row=1, sticky="W")
        self.label_Locus_name = Label(self.label_Gene, text="")
        self.label_Locus_name.grid(column=1, row=1, sticky="W")
        self.label_Gene_id = Label(self.label_Gene, text="Gene ID: ")
        self.label_Gene_id.grid(column=0, row=2, sticky="W")
        self.label_Gene_id_name = Label(self.label_Gene, text="")
        self.label_Gene_id_name.grid(column=1, row=2, sticky="W")
        self.label_Gene_name = Label(self.label_Gene, text="Gene Name: ")
        self.label_Gene_name.grid(column=0, row=3, sticky="W")
        self.label_Gene_name_name = Label(self.label_Gene, text="")
        self.label_Gene_name_name.grid(column=1, row=3, sticky="W")
        self.label_Gene_Size = Label(self.label_Gene, text="Size: ")
        self.label_Gene_Size.grid(column=0, row=4, sticky="W")
        self.label_Gene_Size_name = Label(self.label_Gene, text="")
        self.label_Gene_Size_name.grid(column=1, row=4, sticky="W")
        self.label_Gene_region = Label(self.label_Gene, text="Region: ")
        self.label_Gene_region.grid(column=0, row=5, sticky="W")
        self.label_Gene_region_name = Label(self.label_Gene, text="")
        self.label_Gene_region_name.grid(column=1, row=5, sticky="W")
        self.label_protein_name = Label(self.label_Gene, text="Protein: ")
        self.label_protein_name.grid(column=0, row=6, sticky="W")
        self.label_protein_name_name = Label(self.label_Gene, text="")
        self.label_protein_name_name.grid(column=1, row=6, sticky="W")
        self.label_protein_id = Label(self.label_Gene, text="Protein ID: ")
        self.label_protein_id.grid(column=0, row=7, sticky="W")
        self.label_protein_id_name = Label(self.label_Gene, text="")
        self.label_protein_id_name.grid(column=1, row=7, sticky="W")

        # ===================================
        # Add a PCR Template Frame to Tab 2
        # ===================================
        self.label_frame_0 = LabelFrame(self.tab2, text="PCR Template", fg="red")
        self.label_frame_0.grid(column=0, row=3, sticky='W')

        # ===================================
        # Add a Selection  Frame to Pcr Template Frame
        # ===================================
        self.Selecton_frame = LabelFrame(self.label_frame_0)
        self.Selecton_frame.place(x=240, y=-5, width=398, height=25)
        # ====================================================
        # Add Selection option Menu  to Selection Frame
        # ====================================================
    def Selection_quary(self):
        self.organism_clicked = StringVar()
        self.organism_drop = ttk.Combobox(self.Selecton_frame, width= 15, textvariable= self.organism_clicked)
        self.organism_drop["values"]= self.quary_data()
        self.organism_drop.bind("<<ComboboxSelected>>",self.Species_quary)
        self.organism_drop.set("Genus")
        self.organism_drop.place(x=-3, y=-5, width=100, height=30)

    # ====================================================
    # *******Quary Database ********
    # ====================================================
    def quary_data(self):
        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')

        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()
        self.genus_quary = []
        self.cur.execute("SELECT Genus FROM  Organism")
        self.quary = self.cur.fetchall()
        for geneus in self.quary:
            if geneus[0] not in self.genus_quary:
                self.genus_quary.append(geneus[0])
        return self.genus_quary
        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()
        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()

    # ====================================================
    # ******* Species Quary Database ********
    # ====================================================
    def Species_quary(self, event=None):
        self.Species_clicked = StringVar()
        self.Species_drop = ttk.Combobox(self.Selecton_frame, width=15, textvariable=self.Species_clicked)
        self.Species_drop.set("Species")
        self.Species_drop.place(x=95, y=-5, width=100, height=30)
        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')
        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()
        self.species_quary = []
        self.genus_quary_name = self.organism_drop.get()
        if not self.genus_quary_name == "Genus":
            self.cur.execute("""SELECT Species FROM  Organism WHERE  Genus = ?""", (self.genus_quary_name,))
            self.quary_species = self.cur.fetchall()
            for species in self.quary_species:
                if species[0] not in self.species_quary:
                    self.species_quary.append(species[0])
            print(self.species_quary)
        self.Species_drop["values"] = self.species_quary
        self.Species_drop.bind("<<ComboboxSelected>>", self.Sero_quary)
        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()
        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()

    # ====================================================
    # ******* SeroVAr Quary Database ********
    # ====================================================
    def Sero_quary(self, event=None):
        self.Sero_clicked = StringVar()
        self.Sero_drop = ttk.Combobox(self.Selecton_frame, width=15, textvariable=self.Sero_clicked)
        self.Sero_drop.set("Serovar")
        self.Sero_drop.place(x=195, y=-5, width=100, height=30)
        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')
        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()
        self.serovar_quary = []
        self.species_quary_name = self.Species_drop.get()
        if not self.species_quary_name == "Species":
            self.cur.execute("""SELECT Serovar FROM  Organism WHERE  Genus = ? AND Species = ?""",
                             (self.genus_quary_name, self.species_quary_name))
            self.quary_serovar = self.cur.fetchall()
            for serovar in self.quary_serovar:
                if serovar[0] not in self.serovar_quary:
                    self.serovar_quary.append(serovar[0])
            print(self.serovar_quary)
        self.Sero_drop["values"] = self.serovar_quary
        self.Sero_drop.bind("<<ComboboxSelected>>", self.Gene_quary)
        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()
        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()

    # ====================================================
    # ******* Gene Quary Database ********
    # ===================================================
    def Gene_quary(self, event=None):
        self.Gene_clicked = StringVar()
        self.Gene_drop = ttk.Combobox(self.Selecton_frame, width=15, textvariable=self.Gene_clicked)
        self.Gene_drop.set("Gene")
        self.Gene_drop.place(x=295, y=-5, width=100, height=30)
        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')
        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()
        self.gene_quary = []
        self.Gene_quary_name = self.Sero_drop.get()
        if not self.Gene_quary_name == "Gene":
            self.cur.execute(
                """SELECT Gene_Name FROM  Gene, Organism WHERE  Genus = ? AND Species = ? AND Serovar = ? AND  Organism.Accession = Gene.Accession """,
                (self.genus_quary_name, self.species_quary_name, self.Gene_quary_name,))
            self.quary_gene = self.cur.fetchall()
            for Gene in self.quary_gene:
                if Gene not in self.gene_quary:
                    self.gene_quary.append(Gene[0])
            print(self.gene_quary)
        self.Gene_drop["values"] = self.gene_quary
        self.Gene_drop.bind("<<ComboboxSelected>>", self.Gene_info)
        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()
        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()

    # ====================================================
    # ******* Gene information show ********
    # ====================================================
    # ============================================
    # Create output Gene Information
    # ============================================
    def Gene_info(self, event=None):
        window_info = Tk()
        window_info.title("Gene Information")
        window_info.geometry('700x680+50+25')
        # ================================
        # Disable resizing the output Gene Information
        # ================================
        window_info.resizable(0, 0)

        # =======================================
        # Add "ScrolledText" to output data instance
        # =======================================
        Gene_text = ScrolledText(window_info, width=85, height=40)
        Gene_text.grid(column=0, sticky='WE', columnspan=3)
        # ====================================================
        # *******create Database or connect if found********
        # ====================================================
        conn = sqlite3.connect('vfdb.db')
        # ====================================================
        # *******Dealing with the database********
        # ====================================================
        self.cur = conn.cursor()
        self.gene_quary = []
        self.Gene_info_name = self.Gene_drop.get()
        if not self.Gene_quary_name == "Gene":
            self.cur.execute("""SELECT Sequence  FROM  Gene WHERE  Gene_Name = ?  """,(self.Gene_info_name,))
            self.info_gene = self.cur.fetchall()
            self.scr.delete('0.0', END)
            self.scr.insert(INSERT, self.info_gene[0])
            self.cur.execute("""SELECT *  FROM  Gene, Organism WHERE  Gene_Name = ?  AND  Organism.Accession = Gene.Accession """, (self.Gene_info_name,))
            self.info_gene_1 = self.cur.fetchall()
            print(self.info_gene_1[0])
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "  \n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, " Ogranism  Information  ")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "==============================")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT,"{0:15s}  {1:25s}  {2:9s} {3:2d}".format(" Accession     :", self.info_gene_1[0][7], "Version :", self.info_gene_1[0][18]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT,"{0:15s}  {1:25s}  {2:9s} {3:15s}".format(" Genus         :", self.info_gene_1[0][10], "Species:", self.info_gene_1[0][11]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT,"{0:15s}  {1:25s}  {2:9s} {3:20s}".format(" Sub_Species   :", self.info_gene_1[0][12], "SeroVar :", self.info_gene_1[0][13]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT,"{0:15s}  {1:25s}  {2:9s} {3:20s}".format(" Molecule Types:", self.info_gene_1[0][15], "Taxon  :", self.info_gene_1[0][17]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:46s} ".format(" Strain        :", self.info_gene_1[0][14]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:46s} ".format(" DBLINK        :", self.info_gene_1[0][16][:67]))
        Gene_text.insert(INSERT, "{0:15s}  {1:46s} ".format("               ", self.info_gene_1[0][16][67:]))
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "            **********************************************************")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, " Gene  Information  ")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "==============================")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:30}  {2:9s} {3:15s}".format(" Locus name    :", self.info_gene_1[0][0],"Gene ID :", self.info_gene_1[0][2]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:45s}  ".format(" Gene name     :", self.info_gene_1[0][1]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:30}  {2:9s} {3:<15d}".format(" Region        :", self.info_gene_1[0][4],"Size    :", self.info_gene_1[0][3]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:45s}  ".format(" Protein ID    :", self.info_gene_1[0][6]))
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:15s}  {1:45s}  ".format(" Protein Name  :", self.info_gene_1[0][5] [:66]))
        Gene_text.insert(INSERT, "{0:15s}  {1:45s}  ".format("               ", self.info_gene_1[0][5][66:]))
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "            **********************************************************")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, " Sequence  ")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "==============================")
        Gene_text.insert(INSERT, "  ")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "\n")
        Gene_text.insert(INSERT, "{0:69s} \n ".format(self.info_gene_1[0][8]))
        Gene_text.insert(INSERT, "\n")

        # ====================================================
        # *******the accumulated changes to the database********
        # ====================================================
        conn.commit()
        # ====================================================
        # *******Close the Database********
        # ====================================================
        conn.close()



    # ====================================================
    # ******* Entry Parmeter ********
    # ====================================================

    def Entry_Parmeter(self):
        # ====================================================
        # Add a Entry sequence Textbox  to PCR Template Frame
        # ====================================================
        self.label_Enter_Sequence = Label(self.label_frame_0, text="Enter Sequence or Select Gene:")
        self.label_Enter_Sequence.grid(column=0, row=0, sticky="W", padx=10)
        self.scr = ScrolledText(self.label_frame_0, width=80, height=5)
        self.scr.grid(column=0, sticky='WE', columnspan=3)

        # ===============================================================
        # Add "Selection / simulation  / Blast   Radiobutton" to PCR Template Frame
        # ===============================================================


        # ===============================================================
        # 2- Add "Selection / simulation  / Blast   Radiobutton"
        # ===============================================================

        self.radvar_Sel = IntVar()
        R_Selection_primer = Radiobutton(self.label_frame_0, text=" Primer or Probe Selection", variable=self.radvar_Sel, value=1, command=self.radCall)
        R_Selection_primer.grid(column=0, row=6, sticky="W")
        R_Blast = Radiobutton(self.label_frame_0, text=" Blast", variable=self.radvar_Sel, value=3, command=self.radCall)
        R_Blast.grid(column=2, row=6, sticky="W")
        R_Simulation = Radiobutton(self.label_frame_0, text=" Simulation", variable=self.radvar_Sel, value=2, command=self.radCall)
        R_Simulation.grid(column=1, row=6, sticky="W")

        # ===============================================================
        # Add "Primer / probe Select Checkbutton" to PCR Template Frame
        # ===============================================================
        self.var_p_l = IntVar()
        self.c_left_primer = Checkbutton( self.label_frame_0, text=" Pick left primer", variable=self.var_p_l)
        self.c_left_primer.select()
        self.var_p_r = IntVar()
        self.c_left_primer.grid(column=0, row=7, sticky="W")
        self.c_right_primer = Checkbutton(self.label_frame_0, text=" Pick right primer", variable=self.var_p_r)
        self.c_right_primer.select()
        self.c_right_primer.grid(column=2, row=7, sticky="W")
        self.var_p = IntVar()
        self.c_probe = Checkbutton(self.label_frame_0, text=" Pick hybridization probe", variable=self.var_p, command=self.probe_Call)
        self.c_probe.deselect()
        self.c_probe.grid(column=1, row=7, sticky="W")

        # =========================================
        # Add a "Primer Parameters" Frame to Tab 2
        # =========================================
        label_frame_1 = LabelFrame(self.tab2, text="Primer Parameters", fg="red")
        label_frame_1.grid(column=0, row=4, sticky='W')

        # =========================================================================
        # Adding a "PCR Product Size" Label, and Entry to"Primer Parameters" Frame
        # =========================================================================
        self.label_Product_Size = Label(label_frame_1, text="PCR Product Size")
        self.label_Product_Size.grid(column=0, row=1, sticky="W")
        self.label_Product_Size_1 = Label(label_frame_1, text="Min")
        self.label_Product_Size_1.grid(column=1, row=0)
        self.e_P_1 = StringVar()
        self.entery_Product_Size_1 = Entry(label_frame_1, width=4, textvariable=self.e_P_1)
        self.entery_Product_Size_1.grid(column=1, row=1, sticky='E', padx=100)
        self.e_P_1.set("70")
        self.label_Product_Size_2 = Label(label_frame_1, text="Max")
        self.label_Product_Size_2.grid(column=2, row=0, )
        self.e_P_2 = StringVar()
        self.entery_Product_Size_2 = Entry(label_frame_1, width=4, textvariable=self.e_P_2)
        self.entery_Product_Size_2.grid(column=2, row=1, sticky='E', padx=100)
        self.e_P_2.set("1000")

        # ===============================================================================
        # Adding a "# of primers to return" Label, and Entry to"Primer Parameters" Frame
        # ===============================================================================
        self.label_Primer_return = Label(label_frame_1, text="# of primers to return")
        self.label_Primer_return.grid(column=0, row=2, sticky="W", pady=10)
        self.e_r_1 = StringVar()
        self.entery_Primer_return = Entry(label_frame_1, width=4, textvariable=self.e_r_1)
        self.entery_Primer_return.grid(column=1, row=2, sticky='E', pady=10, padx=100)
        self.e_r_1.set("5")

        # ====================================================================
        # Adding a "Primer Size" Label, and Entry to"Primer Parameters" Frame
        # ====================================================================
        self.label_Primer_Size = Label(label_frame_1, text="Primer Size")
        self.label_Primer_Size.grid(column=0, row=4, sticky="W")
        self.label_Size_1 = Label(label_frame_1, text="Min")
        self.label_Size_1.grid(column=1, row=3)
        self.e_1 = StringVar()
        self.entery_Size_1 = Entry(label_frame_1, width=4, textvariable=self.e_1)
        self.entery_Size_1.grid(column=1, row=4, sticky='E', padx=100)
        self.e_1.set("15")
        self.label_Size_2 = Label(label_frame_1, text="Max")
        self.label_Size_2.grid(column=2, row=3, )
        self.e_2 = StringVar()
        self.entery_Size_2 = Entry(label_frame_1, width=4, textvariable=self.e_2)
        self.entery_Size_2.grid(column=2, row=4, sticky='E', padx=100)
        self.e_2.set("25")
        self.label_Size_3 = Label(label_frame_1, text="Opt")
        self.label_Size_3.grid(column=3, row=3)
        self.e_3 = StringVar()
        self.entery_Size_3 = Entry(label_frame_1, width=4, textvariable=self.e_3)
        self.entery_Size_3.grid(column=3, row=4, sticky='E', padx=10)
        self.e_3.set("20")

        # ===================================================================================
        # Adding a "Primer Melting Temperature" Label, and Entry to"Primer Parameters" Frame
        # ===================================================================================
        self.label_Primer_Temp = Label(label_frame_1, text="Primer Melting Temperature")
        self.label_Primer_Temp.grid(column=0, row=7)
        self.label_Temp_1 = Label(label_frame_1, text="Min")
        self.label_Temp_1.grid(column=1, row=6)
        self.t_1 = StringVar()
        self.entery_Temp_1 = Entry(label_frame_1, width=4, textvariable=self.t_1)
        self.entery_Temp_1.grid(column=1, row=7, sticky='E', padx=100)
        self.t_1.set("57.0")
        self.label_Temp_2 = Label(label_frame_1, text="Max")
        self.label_Temp_2.grid(column=2, row=6, )
        self.t_2 = StringVar()
        self.entery_Temp_2 = Entry(label_frame_1, width=4, textvariable=self.t_2)
        self.entery_Temp_2.grid(column=2, row=7, sticky='E', padx=100)
        self.t_2.set("63.0")
        self.label_Temp_3 = Label(label_frame_1, text="Opt")
        self.label_Temp_3.grid(column=3, row=6)
        self.t_3 = StringVar()
        self.entery_Temp_3 = Entry(label_frame_1, width=4, textvariable=self.t_3)
        self.entery_Temp_3.grid(column=3, row=7, sticky='E', padx=10)
        self.t_3.set("60.0")

        # ===============================================================================
        # Adding a "Primer GC Content (%)" Label, and Entry to"Primer Parameters" Frame
        # ===============================================================================
        self.label_Primer_GC = Label(label_frame_1, text="Primer GC Content (%)")
        self.label_Primer_GC.grid(column=0, row=9, sticky='W')
        self.label_GC_1 = Label(label_frame_1, text="Min")
        self.label_GC_1.grid(column=1, row=8)
        self.gc_1 = StringVar()
        self.entery_GC_1 = Entry(label_frame_1, width=4, textvariable=self.gc_1)
        self.entery_GC_1.grid(column=1, row=9, sticky='E', padx=100)
        self.gc_1.set("20.0")
        self.label_GC_2 = Label(label_frame_1, text="Max")
        self.label_GC_2.grid(column=2, row=8, )
        self.gc_2 = StringVar()
        self.entery_GC_2 = Entry(label_frame_1, width=4, textvariable=self.gc_2)
        self.entery_GC_2.grid(column=2, row=9, sticky='E', padx=100)
        self.gc_2.set("80.0")

        # ===================================================================
        # Adding a "Max Poly_X" Label, and Entry to"Primer Parameters" Frame
        # ===================================================================
        self.label_Poly_X = Label(label_frame_1, text="Max Poly_X")
        self.label_Poly_X.grid(column=0, row=10, sticky='W', pady=10)
        self.poly = StringVar()
        self.entery_Poly_x = Entry(label_frame_1, width=4, textvariable=self.poly)
        self.entery_Poly_x.grid(column=1, row=10, sticky='E', padx=100, pady=10)
        self.poly.set("5")

        # =======================================================================
        # Adding a "Max N Accepted" Label, and Entry to"Primer Parameters" Frame
        # =======================================================================
        self.label_MAx_N = Label(label_frame_1, text="Max N Accepted")
        self.label_MAx_N.grid(column=0, row=11, sticky='W', pady=10)
        self.poly_N = StringVar()
        self.entery_Poly_N = Entry(label_frame_1, width=4, textvariable=self.poly_N)
        self.entery_Poly_N.grid(column=1, row=11, sticky='E', padx=100, pady=10)
        self.poly_N.set("0")

        # =================================================================================
        # Adding a "Max Self Complementarity" Label, and Entry to"Primer Parameters" Frame
        # =================================================================================
        self.label_Max_Self_Compl = Label(label_frame_1, text="Max Self Complementarity")
        self.label_Max_Self_Compl.grid(column=0, row=13, sticky='W')
        self.label_Max_Self_Compl_1 = Label(label_frame_1, text="Any")
        self.label_Max_Self_Compl_1.grid(column=1, row=12)
        self.sc_1 = StringVar()
        self.entery_Max_Self_Compl_1 = Entry(label_frame_1, width=4, textvariable=self.sc_1)
        self.entery_Max_Self_Compl_1.grid(column=1, row=13, sticky='E', padx=100)
        self.sc_1.set("8")
        self.label_Max_Self_Compl_2 = Label(label_frame_1, text="3'")
        self.label_Max_Self_Compl_2.grid(column=2, row=12, )
        self.sc_2 = StringVar()
        self.entery_Max_Self_Compl_2 = Entry(label_frame_1, width=4, textvariable=self.sc_2)
        self.entery_Max_Self_Compl_2.grid(column=2, row=13, sticky='E', padx=100)
        self.sc_2.set("3")

        # =================================================================================
        # Adding a "Max Pair Complementarity" Label, and Entry to"Primer Parameters" Frame
        # =================================================================================
        self.label_Max_Pair_Compl = Label(label_frame_1, text="Max Pair Complementarity")
        self.label_Max_Pair_Compl.grid(column=0, row=15, sticky='W')
        self.label_Max_Pair_Compl_1 = Label(label_frame_1, text="Any")
        self.label_Max_Pair_Compl_1.grid(column=1, row=14)
        self.pc_1 = StringVar()
        self.entery_Max_Pair_Compl_1 = Entry(label_frame_1, width=4, textvariable=self.pc_1)
        self.entery_Max_Pair_Compl_1.grid(column=1, row=15, sticky='E', padx=100)
        self.pc_1.set("8")
        self.label_Max_Pair_Compl_2 = Label(label_frame_1, text="3'")
        self.label_Max_Pair_Compl_2.grid(column=2, row=14, )
        self.pc_2 = StringVar()
        self.entery_Max_Pair_Compl_2 = Entry(label_frame_1, width=4, textvariable=self.pc_2)
        self.entery_Max_Pair_Compl_2.grid(column=2, row=15, sticky='E', padx=100)
        self.pc_2.set("3")

        # =========================================
        # Add a "Buffer Condition" Frame to Tab 2
        # =========================================
        label_frame_2 = LabelFrame(self.tab2, text="Buffer Condition", fg="red")
        label_frame_2.grid(column=0, row=5, sticky='W')

        # ===================================================================================
        # Adding a "Condition" Label, and Entry to"Buffer Condition" Frame
        # ===================================================================================
        self.label_con_mon_cat = Label(label_frame_2, text=" Concentration of Monovalent Cations ")
        self.label_con_mon_cat.grid(column=0, row=0, sticky='W')
        self.c_m_1 = StringVar()
        self.entery_con_mon_cat = Entry(label_frame_2, width=4, textvariable=self.c_m_1)
        self.entery_con_mon_cat.grid(column=1, row=0, sticky='E', padx=55)
        self.c_m_1.set("50.0")

        self.label_con_dNTP = Label(label_frame_2, text=" Concentration of dNTPs ")
        self.label_con_dNTP.grid(column=2, row=0, sticky='W')
        self.c_dNTP_1 = StringVar()
        self.entery_con_dNTP = Entry(label_frame_2, width=4, textvariable=self.c_dNTP_1)
        self.entery_con_dNTP.grid(column=3, row=0, sticky='E', padx=55)
        self.c_dNTP_1.set("0.6")

        self.label_con_di_cat = Label(label_frame_2, text=" Concentration of Divalent Cations ")
        self.label_con_di_cat.grid(column=0, row=1, sticky='W')
        self.c_d_1 = StringVar()
        self.entery_con_di_cat = Entry(label_frame_2, width=4, textvariable=self.c_d_1)
        self.entery_con_di_cat.grid(column=1, row=1, sticky='E', padx=55)
        self.c_d_1.set("1.5")

        self.label_con_AO = Label(label_frame_2, text=" Annealing Oligo Concentration")
        self.label_con_AO.grid(column=2, row=1, sticky='W')
        self.c_AO_1 = StringVar()
        self.entery_con_AO = Entry(label_frame_2, width=4, textvariable=self.c_AO_1)
        self.entery_con_AO.grid(column=3, row=1, sticky='E', padx=55)
        self.c_AO_1.set("50.0")

        # ========================================================
        # Add "Pick Primer / probe button" to PCR Template Frame
        # ========================================================
        button = Button(self.label_frame_0, text="Pick Primer", command=self.second_window)
        button.grid(column=0, row=3, sticky='WE')

        # ========================================================
        # Add "Comparative analysis button" to PCR Template Frame
        # ========================================================
        button = Button(self.label_frame_0, text="Compare", command=self.Compare_window)
        button.grid(column=2, row=3, sticky='WE')

        # ========================================================
        # Add "retrieve  button" to browser Frame
        # ========================================================
        self.label_retrieve_button = Button(self.label_open_file_frame, text="retrieve",state='disabled', command=self.retrive_gene)
        self.label_retrieve_button.grid(column=3, row=0, sticky='WE')


        # ========================================================
        # Add "Database file browser button" to browser Frame
        # ========================================================
        self.label_open_file_button = Button(self.label_open_file_frame, text="Browse A File", state='disabled', command=self.open_file)
        self.label_open_file_button.grid(column=3, row=1, sticky='WE')

        # =========================
        # Create output function
        # =========================
        def radCall(self):
            radSel = self.radvar_Sel.get()
            if radSel == 2 or radSel == 3:
                self.c_left_primer.configure(state='disabled')
                self.c_right_primer.configure(state='disabled')
                self.c_probe.configure(state='disabled')


    def Compare_window(self):
        #input = self.scr.get("1.0", "end-1c")
        # =========================
        # Create output data instance
        # =========================
        compare_window = Tk()
        compare_window.title("Gene Comparative Analysis")
        compare_window.geometry('1300x1000+50+25')

        # ================================
        # Disable resizing the output data instance
        # ================================
        compare_window.resizable(0, 0)

        # =======================================
        # Add "ScrolledText" to output data instance
        # =======================================
        compar_text = ScrolledText(compare_window, width=180, height=40)
        compar_text.grid(column=0, row=0, padx=10)
        #result_handle = open("results.xml")    this line to chick the function without need to the internet to save time
        result_handle = NCBIWWW.qblast("blastn", "nt", self.scr.get("1.0", "end-1c"))
        E_VALUE_THRESH = 0.001
        for record in NCBIXML.parse(result_handle):
            compar_text.insert(INSERT, "                     \n")
            compar_text.insert(INSERT, "                                                         This gene is matching with the below Aligment \n")
            compar_text.insert(INSERT, "                                                                   **************************          \n")
            compar_text.insert(INSERT, "                     \n")
            if record.alignments:
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            compar_text.insert(INSERT, "  match: %s \n" % align.title[:150])
                            break


    def second_window(self):
        input = self.scr.get("1.0", "end-1c")
        # =========================
        # Create output data instance
        # =========================
        window = Tk()
        window.title("Pick Primer")
        window.geometry('1228x680+50+25')

        # ================================
        # Disable resizing the output data instance
        # ================================
        window.resizable(0, 0)

        # =======================================
        # Add "ScrolledText" to output data instance
        # =======================================
        text = ScrolledText(window, width=150, height=40)
        list_cycl_prd = [st]
        if self.radvar_Sel.get() == 1:
            if self.var_p_l.get() == 1 and self.var_p_r.get() == 1:
                primer = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'mohamed',
                        'SEQUENCE_TEMPLATE': input,
                        'SEQUENCE_EXCLUDED_REGION': [0, 0]
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_NUM_RETURN': int(self.e_r_1.get()),
                        'PRIMER_OPT_SIZE': int(self.e_3.get()),
                        'PRIMER_MIN_SIZE': int(self.e_1.get()),
                        'PRIMER_MAX_SIZE': int(self.e_2.get()),
                        'PRIMER_OPT_TM': float(self.t_3.get()),
                        'PRIMER_MIN_TM': float(self.t_1.get()),
                        'PRIMER_MAX_TM': float(self.t_2.get()),
                        'PRIMER_MIN_GC': float(self.gc_1.get()),
                        'PRIMER_MAX_GC': float(self.gc_2.get()),
                        'PRIMER_MAX_POLY_X': int(self.poly.get()),
                        'PRIMER_SALT_MONOVALENT': float(self.c_m_1.get()),
                        'PRIMER_SALT_DIVALENT': float(self.c_d_1.get()),
                        'PRIMER_DNTP_CONC': float(self.c_dNTP_1.get()),
                        'PRIMER_DNA_CONC': float(self.c_AO_1.get()),
                        'PRIMER_MAX_NS_ACCEPTED': int(self.poly_N.get()),
                        'PRIMER_MAX_SELF_ANY': int(self.sc_1.get()),
                        'PRIMER_MAX_SELF_END': int(self.sc_2.get()),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(self.pc_1.get()),
                        'PRIMER_PAIR_MAX_COMPL_END': int(self.pc_2.get()),
                        'PRIMER_PRODUCT_SIZE_RANGE': [int(self.e_P_1.get()), int(self.e_P_2.get())], }))
                for x in range(0, int(self.e_r_1.get())):
                    start_lift = primer["PRIMER_LEFT_" + str(x)][0]
                    len_primer_lift = primer['PRIMER_LEFT_' + str(x)][1]
                    stop_lift = start_lift + len_primer_lift - 1
                    primer_tm_lift = primer["PRIMER_LEFT_" + str(x) + "_TM"]
                    primer_gc_lift = primer["PRIMER_LEFT_" + str(x) + "_GC_PERCENT"]
                    primer_any_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_ANY_TH"]
                    primer_3_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_END_TH"]
                    primer_seq_lift = primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"]
                    start_right = primer["PRIMER_RIGHT_" + str(x)][0]
                    len_primer_right = primer["PRIMER_RIGHT_" + str(x)][1]
                    stop_right = start_right - len_primer_right + 1
                    primer_tm_right = primer["PRIMER_RIGHT_" + str(x) + "_TM"]
                    primer_gc_right = primer["PRIMER_RIGHT_" + str(x) + "_GC_PERCENT"]
                    primer_any_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_ANY_TH"]
                    primer_3_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_END_TH"]
                    primer_seq_right = primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"]
                    text.insert(INSERT, " Primer pair  ")
                    text.insert(INSERT, x + 1)
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "  _______________________________________________________________________________________________________________________________________________")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {1:20s}   {2:18s}  {4:6s}  {3:6s}  {5:6s}  {6:6s}  {7:6s}  {8:20s}    {9:20s}  ".format(
                                    "", "Sequence (5'->3')", "Template strand", "Start", "Length", "Stop", "Tm", "GC%",
                                    "Self complementarity", "Self 3' complementarity"))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                    " Forward primer", "Plus", start_lift, len_primer_lift, stop_lift, primer_tm_lift,
                                    primer_gc_lift, primer_any_TH_lift, primer_3_TH_lift, primer_seq_lift))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                    " Reverse primer", "Minus", start_right, len_primer_right, stop_right,
                                    primer_tm_right,
                                    primer_gc_right, primer_any_TH_right, primer_3_TH_right, primer_seq_right))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
            elif self.var_p_l.get() == 1 and self.var_p_r.get() != 1:
                primer = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'mohamed',
                        'SEQUENCE_TEMPLATE': input,
                        'SEQUENCE_EXCLUDED_REGION': [0, 0]
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 0,
                        'PRIMER_NUM_RETURN': int(self.e_r_1.get()),
                        'PRIMER_OPT_SIZE': int(self.e_3.get()),
                        'PRIMER_MIN_SIZE': int(self.e_1.get()),
                        'PRIMER_MAX_SIZE': int(self.e_2.get()),
                        'PRIMER_OPT_TM': float(self.t_3.get()),
                        'PRIMER_MIN_TM': float(self.t_1.get()),
                        'PRIMER_MAX_TM': float(self.t_2.get()),
                        'PRIMER_MIN_GC': float(self.gc_1.get()),
                        'PRIMER_MAX_GC': float(self.gc_2.get()),
                        'PRIMER_MAX_POLY_X': int(self.poly.get()),
                        'PRIMER_SALT_MONOVALENT': float(self.c_m_1.get()),
                        'PRIMER_SALT_DIVALENT': float(self.c_d_1.get()),
                        'PRIMER_DNTP_CONC': float(self.c_dNTP_1.get()),
                        'PRIMER_DNA_CONC': float(self.c_AO_1.get()),
                        'PRIMER_MAX_NS_ACCEPTED': int(self.poly_N.get()),
                        'PRIMER_MAX_SELF_ANY': int(self.sc_1.get()),
                        'PRIMER_MAX_SELF_END': int(self.sc_2.get()),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(self.pc_1.get()),
                        'PRIMER_PAIR_MAX_COMPL_END': int(self.pc_2.get()),
                        'PRIMER_PRODUCT_SIZE_RANGE': [int(self.e_P_1.get()), int(self.e_P_2.get())], }))
                for x in range(0, int(self.e_r_1.get())):
                    start_lift = primer["PRIMER_LEFT_" + str(x)][0]
                    len_primer_lift = primer['PRIMER_LEFT_' + str(x)][1]
                    stop_lift = start_lift + len_primer_lift - 1
                    primer_tm_lift = primer["PRIMER_LEFT_" + str(x) + "_TM"]
                    primer_gc_lift = primer["PRIMER_LEFT_" + str(x) + "_GC_PERCENT"]
                    primer_any_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_ANY_TH"]
                    primer_3_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_END_TH"]
                    primer_seq_lift = primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"]
                    text.insert(INSERT, " Primer pair  ")
                    text.insert(INSERT, x + 1)
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "  _______________________________________________________________________________________________________________________________________________")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {1:20s}   {2:18s}  {4:6s}  {3:6s}  {5:6s}  {6:6s}  {7:6s}  {8:20s}    {9:20s}  ".format(
                                    "", "Sequence (5'->3')", "Template strand", "Start", "Length", "Stop", "Tm", "GC%",
                                    "Self complementarity", "Self 3' complementarity"))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                    " Forward primer", "Plus", start_lift, len_primer_lift, stop_lift, primer_tm_lift,
                                    primer_gc_lift, primer_any_TH_lift, primer_3_TH_lift, primer_seq_lift))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
            elif self.var_p_l.get() != 1 and self.var_p_r.get() == 1:
                primer = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'mohamed',
                        'SEQUENCE_TEMPLATE': input,
                        'SEQUENCE_EXCLUDED_REGION': [0, 0]
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 0,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_NUM_RETURN': int(self.e_r_1.get()),
                        'PRIMER_OPT_SIZE': int(self.e_3.get()),
                        'PRIMER_MIN_SIZE': int(self.e_1.get()),
                        'PRIMER_MAX_SIZE': int(self.e_2.get()),
                        'PRIMER_OPT_TM': float(self.t_3.get()),
                        'PRIMER_MIN_TM': float(self.t_1.get()),
                        'PRIMER_MAX_TM': float(self.t_2.get()),
                        'PRIMER_MIN_GC': float(self.gc_1.get()),
                        'PRIMER_MAX_GC': float(self.gc_2.get()),
                        'PRIMER_MAX_POLY_X': int(self.poly.get()),
                        'PRIMER_SALT_MONOVALENT': float(self.c_m_1.get()),
                        'PRIMER_SALT_DIVALENT': float(self.c_d_1.get()),
                        'PRIMER_DNTP_CONC': float(self.c_dNTP_1.get()),
                        'PRIMER_DNA_CONC': float(self.c_AO_1.get()),
                        'PRIMER_MAX_NS_ACCEPTED': int(self.poly_N.get()),
                        'PRIMER_MAX_SELF_ANY': int(self.sc_1.get()),
                        'PRIMER_MAX_SELF_END': int(self.sc_2.get()),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(self.pc_1.get()),
                        'PRIMER_PAIR_MAX_COMPL_END': int(self.pc_2.get()),
                        'PRIMER_PRODUCT_SIZE_RANGE': [int(self.e_P_1.get()), int(self.e_P_2.get())], }))
                for x in range(0, int(self.e_r_1.get())):
                    start_right = primer["PRIMER_RIGHT_" + str(x)][0]
                    len_primer_right = primer["PRIMER_RIGHT_" + str(x)][1]
                    stop_right = start_right - len_primer_right + 1
                    primer_tm_right = primer["PRIMER_RIGHT_" + str(x) + "_TM"]
                    primer_gc_right = primer["PRIMER_RIGHT_" + str(x) + "_GC_PERCENT"]
                    primer_any_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_ANY_TH"]
                    primer_3_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_END_TH"]
                    primer_seq_right = primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"]
                    text.insert(INSERT, " Primer pair  ")
                    text.insert(INSERT, x + 1)
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "  _______________________________________________________________________________________________________________________________________________")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {1:20s}   {2:18s}  {4:6s}  {3:6s}  {5:6s}  {6:6s}  {7:6s}  {8:20s}    {9:20s}  ".format(
                                    "", "Sequence (5'->3')", "Template strand", "Start", "Length", "Stop", "Tm", "GC%",
                                    "Self complementarity", "Self 3' complementarity"))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                    " Reverse primer", "Minus", start_right, len_primer_right, stop_right,
                                    primer_tm_right,
                                    primer_gc_right, primer_any_TH_right, primer_3_TH_right, primer_seq_right))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
            elif  self.var_p.get() == 1:
                primer = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'mohamed',
                        'SEQUENCE_TEMPLATE': input,
                        'SEQUENCE_EXCLUDED_REGION': [0, 0]
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 1,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_NUM_RETURN': int(self.e_r_1.get()),
                        'PRIMER_OPT_SIZE': int(self.e_3.get()),
                        'PRIMER_MIN_SIZE': int(self.e_1.get()),
                        'PRIMER_MAX_SIZE': int(self.e_2.get()),
                        'PRIMER_OPT_TM': float(self.t_3.get()),
                        'PRIMER_MIN_TM': float(self.t_1.get()),
                        'PRIMER_MAX_TM': float(self.t_2.get()),
                        'PRIMER_MIN_GC': float(self.gc_1.get()),
                        'PRIMER_MAX_GC': float(self.gc_2.get()),
                        'PRIMER_MAX_POLY_X': int(self.poly.get()),
                        'PRIMER_SALT_MONOVALENT': float(self.c_m_1.get()),
                        'PRIMER_SALT_DIVALENT': float(self.c_d_1.get()),
                        'PRIMER_DNTP_CONC': float(self.c_dNTP_1.get()),
                        'PRIMER_DNA_CONC': float(self.c_AO_1.get()),
                        'PRIMER_MAX_NS_ACCEPTED': int(self.poly_N.get()),
                        'PRIMER_MAX_SELF_ANY': int(self.sc_1.get()),
                        'PRIMER_MAX_SELF_END': int(self.sc_2.get()),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(self.pc_1.get()),
                        'PRIMER_PAIR_MAX_COMPL_END': int(self.pc_2.get()),
                        'PRIMER_PRODUCT_SIZE_RANGE': [int(self.e_P_1.get()), int(self.e_P_2.get())], }))
                for x in range(0, int(self.e_r_1.get())):
                    start_lift = primer["PRIMER_LEFT_" + str(x)][0]
                    primer_seq_lift = primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"]
                    start_right = primer["PRIMER_RIGHT_" + str(x)][0]
                    primer_seq_right_org = primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"]
                    primer_seq_right= str(Seq(primer_seq_right_org).reverse_complement())  #to gave the end of probe we convert the str to seq to reverse complemet then to str
                    primer_seq_internal = primer["PRIMER_INTERNAL_" + str(x) + "_SEQUENCE"]
                    probe_lenth = start_right - start_lift
                    probe_seq = input[start_lift:start_right+1]
                    text.insert(INSERT, " Probe  ")
                    text.insert(INSERT, x + 1)
                    text.insert(INSERT, "\n")
                    text.insert(INSERT,
                                "  _______________________________________________________________________________________________________________________________________________")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "{0:14s}  {1:5s}   {2:25s}  {3:5s}  {4:25s}  {5:5s}  {6:25s}  {7:5s}  {8:14s}".format("Probe Starting", " ", "Probe Start Sequence", " ", "Probe Internal Sequence", " ", "Probe End Sequence", " ", "Probe Ending"))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "{0:^14d}  {1:5s}   {2:25s}  {3:5s}  {4:25s}  {5:5s}  {6:25s}  {7:5s}  {8:^14d}".format(start_lift, " ", primer_seq_lift, " ", primer_seq_internal, " ", primer_seq_right, " ", start_right))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "{0:25s}  {1:6d}  ".format("Probe Lenth   :", probe_lenth))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "{0:25s}  ".format("Probe Sequence:"))
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, probe_seq)
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "                                             *******************************************")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
                    text.insert(INSERT, "  ")
                    text.insert(INSERT, "\n")
        elif self.radvar_Sel.get() == 2:
            primer = (primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'mohamed',
                    'SEQUENCE_TEMPLATE': input,
                    'SEQUENCE_EXCLUDED_REGION': [0, 0]
                },
                {
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_INTERNAL_OLIGO': 0,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': int(self.e_r_1.get()),
                        'PRIMER_OPT_SIZE': int(self.e_3.get()),
                        'PRIMER_MIN_SIZE': int(self.e_1.get()),
                        'PRIMER_MAX_SIZE': int(self.e_2.get()),
                        'PRIMER_OPT_TM': float(self.t_3.get()),
                        'PRIMER_MIN_TM': float(self.t_1.get()),
                        'PRIMER_MAX_TM': float(self.t_2.get()),
                        'PRIMER_MIN_GC': float(self.gc_1.get()),
                        'PRIMER_MAX_GC': float(self.gc_2.get()),
                        'PRIMER_MAX_POLY_X': int(self.poly.get()),
                        'PRIMER_SALT_MONOVALENT': float(self.c_m_1.get()),
                        'PRIMER_SALT_DIVALENT': float(self.c_d_1.get()),
                        'PRIMER_DNTP_CONC': float(self.c_dNTP_1.get()),
                        'PRIMER_DNA_CONC': float(self.c_AO_1.get()),
                        'PRIMER_MAX_NS_ACCEPTED': int(self.poly_N.get()),
                        'PRIMER_MAX_SELF_ANY': int(self.sc_1.get()),
                        'PRIMER_MAX_SELF_END': int(self.sc_2.get()),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(self.pc_1.get()),
                        'PRIMER_PAIR_MAX_COMPL_END': int(self.pc_2.get()),
                        'PRIMER_PRODUCT_SIZE_RANGE': [int(self.e_P_1.get()), int(self.e_P_2.get())], }))
            for x in range(0, int(self.e_r_1.get())):
                start_lift = primer["PRIMER_LEFT_" + str(x)][0]
                len_primer_lift = primer['PRIMER_LEFT_' + str(x)][1]
                stop_lift = start_lift + len_primer_lift - 1
                primer_tm_lift = primer["PRIMER_LEFT_" + str(x) + "_TM"]
                primer_gc_lift = primer["PRIMER_LEFT_" + str(x) + "_GC_PERCENT"]
                primer_any_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_ANY_TH"]
                primer_3_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_END_TH"]
                primer_seq_lift = primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"]
                start_right = primer["PRIMER_RIGHT_" + str(x)][0]
                len_primer_right = primer["PRIMER_RIGHT_" + str(x)][1]
                stop_right = start_right - len_primer_right + 1
                primer_tm_right = primer["PRIMER_RIGHT_" + str(x) + "_TM"]
                primer_gc_right = primer["PRIMER_RIGHT_" + str(x) + "_GC_PERCENT"]
                primer_any_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_ANY_TH"]
                primer_3_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_END_TH"]
                primer_seq_right = primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"]
                template = Dseqrecord(input)
                p1 = Seq(primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"])
                p2 = Seq(primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"])
                cycl_prd = pcr(p1, p2, template)
                list_cycl_prd.append([cycl_prd])
                text.insert(INSERT, " Primer pair  ")
                text.insert(INSERT, x + 1)
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "  _______________________________________________________________________________________________________________________________________________")
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "{0:15s}  {1:20s}   {2:18s}  {4:6s}  {3:6s}  {5:6s}  {6:6s}  {7:6s}  {8:20s}    {9:20s}  ".format(
                                "", "Sequence (5'->3')", "Template strand", "Start", "Length", "Stop", "Tm", "GC%",
                                "Self complementarity", "Self 3' complementarity"))
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                " Forward primer", "Plus", start_lift, len_primer_lift, stop_lift, primer_tm_lift,
                                primer_gc_lift, primer_any_TH_lift, primer_3_TH_lift, primer_seq_lift))
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                " Reverse primer", "Minus", start_right, len_primer_right, stop_right, primer_tm_right,
                                primer_gc_right, primer_any_TH_right, primer_3_TH_right, primer_seq_right))
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "                        ***************************          ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, cycl_prd.figure())
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, cycl_prd.program())
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "************************************************************************************")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
        elif self.radvar_Sel.get() == 3:
            primer = (primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'mohamed',
                    'SEQUENCE_TEMPLATE': input,
                    'SEQUENCE_EXCLUDED_REGION': [0, 0]
                },
                {
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_INTERNAL_OLIGO': 0,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_NUM_RETURN': int(self.e_r_1.get()),
                        'PRIMER_OPT_SIZE': int(self.e_3.get()),
                        'PRIMER_MIN_SIZE': int(self.e_1.get()),
                        'PRIMER_MAX_SIZE': int(self.e_2.get()),
                        'PRIMER_OPT_TM': float(self.t_3.get()),
                        'PRIMER_MIN_TM': float(self.t_1.get()),
                        'PRIMER_MAX_TM': float(self.t_2.get()),
                        'PRIMER_MIN_GC': float(self.gc_1.get()),
                        'PRIMER_MAX_GC': float(self.gc_2.get()),
                        'PRIMER_MAX_POLY_X': int(self.poly.get()),
                        'PRIMER_SALT_MONOVALENT': float(self.c_m_1.get()),
                        'PRIMER_SALT_DIVALENT': float(self.c_d_1.get()),
                        'PRIMER_DNTP_CONC': float(self.c_dNTP_1.get()),
                        'PRIMER_DNA_CONC': float(self.c_AO_1.get()),
                        'PRIMER_MAX_NS_ACCEPTED': int(self.poly_N.get()),
                        'PRIMER_MAX_SELF_ANY': int(self.sc_1.get()),
                        'PRIMER_MAX_SELF_END': int(self.sc_2.get()),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(self.pc_1.get()),
                        'PRIMER_PAIR_MAX_COMPL_END': int(self.pc_2.get()),
                        'PRIMER_PRODUCT_SIZE_RANGE': [int(self.e_P_1.get()), int(self.e_P_2.get())], }))
            for x in range(0, int(self.e_r_1.get())):
                start_lift = primer["PRIMER_LEFT_" + str(x)][0]
                len_primer_lift = primer['PRIMER_LEFT_' + str(x)][1]
                stop_lift = start_lift + len_primer_lift - 1
                primer_tm_lift = primer["PRIMER_LEFT_" + str(x) + "_TM"]
                primer_gc_lift = primer["PRIMER_LEFT_" + str(x) + "_GC_PERCENT"]
                primer_any_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_ANY_TH"]
                primer_3_TH_lift = primer["PRIMER_LEFT_" + str(x) + "_SELF_END_TH"]
                primer_seq_lift = primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"]
                start_right = primer["PRIMER_RIGHT_" + str(x)][0]
                len_primer_right = primer["PRIMER_RIGHT_" + str(x)][1]
                stop_right = start_right - len_primer_right + 1
                primer_tm_right = primer["PRIMER_RIGHT_" + str(x) + "_TM"]
                primer_gc_right = primer["PRIMER_RIGHT_" + str(x) + "_GC_PERCENT"]
                primer_any_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_ANY_TH"]
                primer_3_TH_right = primer["PRIMER_RIGHT_" + str(x) + "_SELF_END_TH"]
                primer_seq_right = primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"]
                template = Dseqrecord(input)
                p1 = Seq(primer["PRIMER_LEFT_" + str(x) + "_SEQUENCE"])
                p2 = Seq(primer["PRIMER_RIGHT_" + str(x) + "_SEQUENCE"])
                cycl_prd = pcr(p1, p2, template)
                list_cycl_prd.append([cycl_prd])
                Blast = Primer_BLAST(primer_seq_lift, self.genus_quary_name + " " + self.species_quary_name)
                text.insert(INSERT, "\n")
                text.insert(INSERT, " ")
                text.insert(INSERT, " Primer pair  ")
                text.insert(INSERT, x + 1)
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "  _______________________________________________________________________________________________________________________________________________")
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "{0:15s}  {1:20s}   {2:18s}  {4:6s}  {3:6s}  {5:6s}  {6:6s}  {7:6s}  {8:20s}    {9:20s}  ".format(
                                "", "Sequence (5'->3')", "Template strand", "Start", "Length", "Stop", "Tm", "GC%",
                                "Self complementarity", "Self 3' complementarity"))
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                " Forward primer", "Plus", start_lift, len_primer_lift, stop_lift, primer_tm_lift,
                                primer_gc_lift, primer_any_TH_lift, primer_3_TH_lift, primer_seq_lift))
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "{0:15s}  {9:20s}   {1:18s}  {3:<6d}  {2:<6d}  {4:<6d}  {5:<6.2f}  {6:<6.2f}  {7:<20.2f}    {8:<20.2f}      ".format(
                                " Reverse primer", "Minus", start_right, len_primer_right, stop_right, primer_tm_right,
                                primer_gc_right, primer_any_TH_right, primer_3_TH_right, primer_seq_right))
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "                        ***************************          ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, cycl_prd.figure())
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, cycl_prd.program())
                text.insert(INSERT, "\n")
                text.insert(INSERT, "  ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "                        ***************************          ")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "    Blast:")
                text.insert(INSERT, Blast)
                text.insert(INSERT, "\n")
                text.insert(INSERT, "\n")
                text.insert(INSERT, "\n")
                text.insert(INSERT,
                            "************************************************************************************")

        text.config(state="disabled")
        text.grid(column=0, row=0, padx=10)
        if self.radvar_Sel.get() == 2:
            Gel1(list_cycl_prd, gel_len=16).run(plot=True)




if __name__ == '__main__':
    parent = Parent()
    parent.mainloop()