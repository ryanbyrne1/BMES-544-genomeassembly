import sys,os; sys.path.append(os.environ['BMESAHMETDIR']); import bmes
bmes.pipinstall('simplesam')
import simplesam

def BasicMapping(fastqfile, OlsenallaUli, Segrotundu, EcoliK12):
    
    BWAEXE=bmes.bwaexe();
    samfile=fastqfile + '.sam';
    if not bmes.isfileandnotempty(samfile):
        cmd=BWAEXE +' mem "' + OlsenallaUli + '"' ' "' + fastqfile + '"';
        bmes.system_redirecttofile(cmd,samfile);
    print(samfile)