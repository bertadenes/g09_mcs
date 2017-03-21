#!/usr/bin/env python3
import argparse, sys, os, logging
from subprocess import Popen
from time import sleep
import cclib

command = "rung09_quiet"# give your submission script

def tail(f,BLOCK_SIZE = 256):
    f.seek(0,2)
    size = f.tell()
    if size < BLOCK_SIZE: BLOCK_SIZE = size
    f.seek(-BLOCK_SIZE, 2)
    a = f.readlines()
    rv = ""
    for b in a: rv += str(b)
    return rv

def seekOut(inp):
    filename = inp.split('.')[0]
    if os.path.isfile(filename+".out"):
        return True
    else:
        return False

def checkOut(fout):
    with open(fout, 'rb') as outfile:
        if "termination" in tail(outfile):
            if "Error" in tail(outfile):
                logging.info(fout+" did not converged. Moving on.")
                return True
                # logging.info("Gaussian termination in error. Exiting.")
                # sys.exit(0)
            else:
                return True
        else:
            return False

def check_submit_wait(command,filename):
    RUNNING = DONE = False
    if os.path.isfile(filename+".out"):
        RUNNING = True
        DONE = checkOut(filename+".out")
    if not DONE:
        Popen([command,filename+".gjf"])
    while not RUNNING:
        RUNNING = seekOut(filename+".gjf")
        sleep(15)
    while not DONE:
        DONE = checkOut(filename+".out")
        sleep(15)
    logging.info(filename+" is done.")
    return


def processg09output(filePath):
    """Creats a cclib parsed Object with data read from g09 output.
    Args:
        filePath (str): The path to the file to process.
    Returns:
        Object: Parsed cclib object.
        Also extended with level of theory: data.method and data.basis store the last used level
        of theory.
        Check *http://cclib.github.io/data.html* for further info.
    """
    try:
        try:
            output = cclib.parser.ccopen(filePath)
            data = output.parse()
            with open(filePath,'r') as f:
                for line in f:
                    if "Standard basis" in line: data.basis = line.split()[2]
                    if "SCF Done" in line: data.method = line.split(')')[0].split('(')[-1]
                    if "%nproc" in line: data.proc = line.strip()
                    if "%mem" in line: data.mem = line.strip()
                    if "%chk" in line: data.chk = line.strip()
        except AttributeError:
            logging.info(filePath+" is not a g09 output, or not complete.")
            sys.exit(0)
    except IOError:
        logging.info("File cannot be opened."+filePath)
        sys.exit(0)
    return data
    
elements = {1:'H',2:'He',3:'Li',4:'Be',5:'B',6:'C',7:'N',8:'O',9:'F',10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',19:'K',20:'Ca',21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',53:'I',54:'Xe',55:'Cs',56:'Ba',57:'La',58:'Ce',59:'Pr',60:'Nd',61:'Pm',62:'Sm',63:'Eu',64:'Gd',65:'Tb',66:'Dy',67:'Ho',68:'Er',69:'Tm',70:'Yb',71:'Lu',72:'Hf',73:'Ta',74:'W',75:'Re',76:'Os',77:'Ir',78:'Pt',79:'Au',80:'Hg',81:'Tl',82:'Pb',83:'Bi',84:'Po',85:'At',86:'Rn',87:'Fr',88:'Ra',89:'Ac',90:'Th',91:'Pa',92:'U',93:'Np',94:'Pu',95:'Am',96:'Cm',97:'Bk',98:'Cf',99:'Es',100:'Fm',101:'Md',102:'No',103:'Lr',104:'Rf',105:'Db',106:'Sg',107:'Bh',108:'Hs',109:'Mt',110:'Ds',111:'Rg',112:'Cn'}

def createInput(data,filename,coordinates,increments,move):
    com = "#P "+data.method+'/'+data.basis+"int(grid=ultrafine) opt(modred,loose)"
    with open(filename,"w") as f:
        f.write(data.chk+"\n")
        f.write(data.proc+"\n")
        f.write(data.mem+"\n\n")
        f.write(com+"\n\n"+filename+"\n\n")
        f.write(str(data.charge)+' '+str(data.mult)+'\n')
        for j in range(len(data.atomnos)):
            f.write("%s   %f   %f   %f\n" % (elements[data.atomnos[j]],data.atomcoords[-1][j][0],data.atomcoords[-1][j][1],data.atomcoords[-1][j][2]))
        f.write("\n")
        for i in range(len(coordinates)):
            if i == move:
                f.write(coordinates[i]+"S 1 "+str(float(increments[i]))+"\n")
            else:
                f.write(coordinates[i]+"F\n")
        f.write("\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile')
    args = parser.parse_args()
    logging.basicConfig(filename = "logfile", level = logging.INFO, format='%(asctime)s %(message)s')
    MCsection = False
    STDIN = True
    modred = []
    stdinput = []
    with open(args.inputfile,"r") as inp:
        for line in inp:
            if line.strip() == "":
                MCsection = False
            if MCsection:
                modred.append(line.strip())
            if "multicoord" == line.strip():
                STDIN = False
            if STDIN:
                stdinput.append(line)
            if "multicoord" == line.strip():
                MCsection = True


    print(stdinput)
    print(modred)
    coordinates = []
    steps = 0
    increments = []
    structures = []
    for mod in modred:
        if mod.split()[0] == "X":
            pass
        elif mod.split()[0] == "B":
            if len(mod.split()) != 6:
                print("Cannot interpret input:",mod)
                sys.exit(0)
            try:
                steps_n = int(mod.split()[4])
            except:
                print("Cannot interpret input:",mod)
                sys.exit(0)
            if steps == 0:
                steps = steps_n
            elif steps != steps_n:
                print("Coordinates must have identical number of steps.")
                sys.exit(0)
            coord = "B "
            try:
                coord += str(int(mod.split()[1])) + " "
                coord += str(int(mod.split()[2])) + " "
                coordinates.append(coord)
            except:
                print("Cannot interpret input:",mod)
                sys.exit(0)
            try:
                increments.append(float(mod.split()[5]))
            except:
                print("Cannot interpret input:",mod)
                sys.exit(0)    
        elif mod.split()[0] == "A":
            pass
        elif mod.split()[0] == "D":
            pass
        elif mod.split()[0] == "L":
            pass
        else:
            print("Cannot interpret input:",mod)
            sys.exit(0)

    if not os.path.exists("scan"): os.makedirs("scan")
    os.chdir("scan")
    step = 0
    with open(args.inputfile.split('.')[0]+"0.gjf",'w') as f:
        for line in stdinput:
            f.write(line)
        for coord in coordinates:
            f.write(coord+"F\n")
        f.write("\n")
    if os.fork():
        sys.exit(0)

    logging.info("PID: "+str(os.getpid()))
    check_submit_wait(command,args.inputfile.split('.')[0]+"0")

    structures.append(processg09output(args.inputfile.split('.')[0]+"0.out"))
    step += 1
    while step <= steps:
        for i in range(len(coordinates)):
            actfile = args.inputfile.split('.')[0]+str(step)+"_"+str(i+1)
            createInput(structures[-1],actfile+".gjf",coordinates,increments,move=i)
            check_submit_wait(command,actfile)
            if i == 0:
                structures.append(processg09output(actfile+".out"))
            else:
                structures[-1] = processg09output(actfile+".out")
        step += 1
    with open(args.inputfile.split('.')[0]+"_multiscan.xyz","w") as f:
        zeroenergy = structures[0].scfenergies[-1]
        for i in range(len(structures)):
            f.write(str(len(structures[i].atomnos))+"\n")
            f.write("step "+str(i)+" "+str((structures[i].scfenergies[-1]-zeroenergy)*23.06035)+"\n")
            for j in range(len(structures[i].atomnos)):
                f.write("%s   %f   %f   %f\n" % (elements[structures[i].atomnos[j]],structures[i].atomcoords[-1][j][0],structures[i].atomcoords[-1][j][1],structures[i].atomcoords[-1][j][2]))



if __name__ == "__main__":
    main()    
