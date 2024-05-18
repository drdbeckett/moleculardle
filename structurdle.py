import streamlit as st
from streamlit_ketcher import st_ketcher
from datetime import datetime, timezone
from io import BytesIO
import pandas as pd
import random
import base64
import time
import pytz
import csv
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)


      #####################
      #     A DAB app     #
      # Started 5/10/2024 #
      #####################

state = st.session_state

# TODO: Change the colors automatically, maybe dark and purple?

# initialize some session state values ####################################
# just get the EST date
if 'today' not in state:
    local_tz = "US/Eastern"
    state.today = pytz.timezone(local_tz).localize(datetime.today())

if 'guessnum' not in state:
    state.guessnum = 0

# Block for Win/Loss condition
if 'LockOut' not in state:
    state.LockOut = False
if 'Won' not in state:
    state.Won = False
if 'Lost' not in state:
    state.Lost = False

# Endless mode management, "FirstEndless" is for making sure the sketcher guess doesn't populate
if 'Endless' not in state:
    state.Endless = False
if 'NewEndless' not in state:
    state.NewEndless = False
if 'FirstEndless' not in state:
    state.FirstEndless=False
# molecular weight for daily guess (never change, don't you ever change it)
if 'DailyMW' not in state:
    state.DailyMW = 500.0
# molecular weight for endless, should change with a slider
if 'EndlessMW' not in state:
    state.EndlessMW = 500.0


# store list of guessed SMILES strings to quickly check if a guess is unique
if 'guesses' not in state:
    state.guesses = []
if 'outdf' not in state:
    state.outdf = pd.DataFrame({"Guess Number": [],
                                           "Tanimoto": [],
                                           "MCS": []})

st.set_page_config(layout='wide')
#st.title("STRUCTURDLE or maybe DRUGDLE or maybe MOLECULARDLE")
st.title("STRUCTURDLE or maybe MOLECULARDLE?")
if not state.Endless:
    dailytitle="Daily Puzzle: "+str(state.today.month)+"/"+str(state.today.day)+"/"+str(state.today.year)
    st.title(dailytitle)
if state.Endless:
    st.title("🚨 ENDLESS MODE ACTIVATED 🚨")
# TODO: Implement a molecular weight slider for endless mode

col1, col2 = st.columns([2,1])


# Initial functions and target line state variable setting ####
def read_csv(filename):
    lines = []
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=',', quotechar='"')
        for row in reader:
            lines.append(row)
    return lines

# filter compounds by molecular weight and if the SMILES is valid
def filter_compounds(lines, MW):
    filtered_data = []
    for row in lines:
        smiles = row[3]  # Assuming the SMILES string is in the fourth column
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None and Chem.rdMolDescriptors.CalcExactMolWt(mol) <= MW:
            filtered_data.append(row)
    return filtered_data

# get target line from the file subject to molecular weight cutoff and random number seed
def target_acquisition(MW,rseed):
    lines = read_csv('dsmiles_cleaned.csv')
    MW_filtered = filter_compounds(lines, MW)
    random.seed(rseed)
    targetindex=random.randrange(1,len(MW_filtered))
    print("Weight: ",MW)
    print("Length of list: ",len(MW_filtered))
    return MW_filtered[targetindex-1]

# Generate the line number for grabbing the target from the current date     
# and pull out the full line as an array to parse at our leisure.
# 0 = DBID, 1 = name, 3 = SMILES, 9 = summary
if 'targetline' not in state:
    #dateseed=str(date.today().day) + str(date.today().year) + str(date.today().month)
    dateseed=str(state.today.day) + str(state.today.year) + str(state.today.month)
    state.targetline = target_acquisition(state.DailyMW,dateseed)
    # debug block 
    #print(state.targetline[1])
    #print(state.targetline[3])
    #print(state.targetline[9])
    #quit()
####################################################################

####### Function definitions ######################################
# For visualizing the MCS
def view_mcs(targetm,guessm):
    rgba_color = (0.0, 0.0, 1.0, 0.2) # transparent blue

    # Actual MCS call
    # the default is very loosey-goosey
    #mcs = rdFMCS.FindMCS([targetm,guessm])
    # I think this one is just right
    mcs = rdFMCS.FindMCS([targetm,guessm], matchValences=True, ringMatchesRingOnly=True)
    # This is too strict - I don't want people to not think something's there if it's actually a ring and they drew a methyl
    #mcs = rdFMCS.FindMCS([targetm,guessm], matchValences=True, ringMatchesRingOnly=True, completeRingsOnly=True)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    match1 = guessm.GetSubstructMatch(mcs_mol)

    # Getting the bonds to highlight
    bonds = []
    for bond in mcs_mol.GetBonds():
       aid1 = match1[bond.GetBeginAtomIdx()]
       aid2 = match1[bond.GetEndAtomIdx()]
       bonds.append(guessm.GetBondBetweenAtoms(aid1,aid2).GetIdx())
    mcs_pil = Draw.MolToImage(guessm, highlightAtoms=match1, highlightBonds=bonds, highlightColor=rgba_color)
    #mcs_pil = Draw.MolToImage(mcs_mol)


    # convert the png pil to a data url of a jpeg
    buffered = BytesIO()
    mcs_pil.save(buffered, format="JPEG")
    mcs_b64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return 'data:image/jpeg;base64,' + mcs_b64

# Printing the emoji string for winners/losers
def emojify():
    emojistring=""
    for i in state.outdf['Tanimoto']:
        if float(i) < 0.3:
            emojistring = emojistring+"🟥"
        elif float(i) < 0.6:
            emojistring = emojistring+"🟨"
        elif float(i) <= 1.0:
            emojistring = emojistring+"🟩"
    if state.Won:
        emojistring = emojistring+"🧪"
    if state.Lost:
        emojistring = emojistring+"☠️"
    return emojistring

# for clearing variables in Endless Mode
def clean_slate():
    state.Won = False
    state.Lost = False
    state.LockOut = False
    state.Endless = True
    state.NewEndless = True
    state.FirstEndless = True
    state.guessnum = 0
    state.outdf = pd.DataFrame({"Guess Number": [],
                                "Tanimoto": [],
                                "MCS": []})

#######################################################
# Target definition and initialization
### debug targets
#target = "NC1CCCCC1"
#target = "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4"
#target = "O=C(N(C)C1=C2N(C)C=N1)N(C)C2=O"
#target = "C[N+]1(C)[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)C1=CC=CC=C1"

# Get the target SMILES string from the csv file
if state.NewEndless:      
    timeseed=int(time.time()*1000)
    state.targetline = target_acquisition(state.EndlessMW,timeseed)
    state.NewEndless=False
    state.FirstEndless=True
target=state.targetline[3].strip('\"')
targetm = Chem.MolFromSmiles(target)
targetHAC = rdMolDescriptors.CalcNumHeavyAtoms(targetm)
targetNumHD = rdMolDescriptors.CalcNumHBD(targetm)
targetNumHA = rdMolDescriptors.CalcNumHBA(targetm)
#targetRC = rdMolDescriptors.CalcNumRings(targetm)
targetArRC = rdMolDescriptors.CalcNumAromaticRings(targetm)
targetAlRC = rdMolDescriptors.CalcNumAliphaticRings(targetm)
targetfp = FingerprintMols.FingerprintMol(targetm)
targetformula = rdMolDescriptors.CalcMolFormula(targetm)


guess = ""
validguess = False

# for the winners or losers
if state.LockOut:
    # TODO: Need to define a function that changes the cookie depending on the score for the day, it should store the state and the day
    # then if the date does not match it pushes the old state to a counter so we can keep track of averages. Ideally we'd keep track of 
    # just a few stats:
              # 1. Today's result (initiate lockout if the game has already been played today)
              # 2. Today's working score for the daily puzzle (do this LAST, might be heavy)
              # 3. Win streak (maybe a lose streak and a daily streak later, keep it simple for the moment)
              # 4. Average number of guesses per win on Daily
              # 5. Average number of wins per Daily play
    if state.Won:
        st.write("You got it on guess ",str(state.guessnum),"! The answer is ",state.targetline[1].strip('\"'))
        link='https://go.drugbank.com/drugs/'+state.targetline[0].strip('\"')
        st.write(link)
        st.write(state.targetline[9].strip('\"'))
        st.image(view_mcs(state.FinalGuessm,targetm))
        if not state.Endless:  
            st.write("Copy the emoji string to show off to your friends, colleagues, and enemies!")
            st.write("Structurdle ",str(state.today.month),"/",str(state.today.day),"/",str(state.today.year),": ", emojify())

    if state.Lost:
        # TODO: Move the image up and try to make it larger  
        st.write("Better luck next time! Here's how close you got.")
        st.write("The answer is ",state.targetline[1].strip('\"'))
        link='https://go.drugbank.com/drugs/'+state.targetline[0].strip('\"')
        st.write(link)
        st.write(state.targetline[9].strip('\"'))
        st.image(view_mcs(state.FinalGuessm,targetm))
        if not state.Endless:  
            st.write("Copy the emoji string to demonstrate how hard you tried before tapping out!")
            st.write("Structurdle ",str(state.today.month),"/",str(state.today.day),"/",str(state.today.year),": ", emojify())

    # reset button that initiates endless mode
    if st.button(":green-background[♾️ Continue in Endless Mode? ♾️]", type="secondary"):
        clean_slate()
        st.rerun()

    state.EndlessMW = st.slider("Endless Mode Molecular Weight Cutoff", 50, 1000, int(state.EndlessMW))
###########################

# Get input structure and properties
# TODO: style these writes
if not state.LockOut:
    with col1:
        st.write("Guess the drug (or drug-like compound)!")
        st.write("Target empirical formula:", targetformula)
        guess = st_ketcher()

# get properties from the input guess
if guess and not state.FirstEndless: 
    if guess not in state.guesses:
        # TODO: add bit that checks against molecular formula
        state.guesses.append(guess)
        state.guessnum = state.guessnum+1
        validguess = True
    else:
        st.write("You already tried that one! Try something else.")
    guessm = Chem.MolFromSmiles(guess)
    guessfp = FingerprintMols.FingerprintMol(guessm)
    guessHAC = rdMolDescriptors.CalcNumHeavyAtoms(guessm)
    guessNumHD = rdMolDescriptors.CalcNumHBD(guessm)
    guessNumHA = rdMolDescriptors.CalcNumHBA(guessm)
    #guessRC = rdMolDescriptors.CalcNumRings(guessm)
    guessArRC = rdMolDescriptors.CalcNumAromaticRings(guessm)
    guessAlRC = rdMolDescriptors.CalcNumAliphaticRings(guessm)

# The give up button
if not state.LockOut:
    with col2:
        if st.button("☠️Give up?☠️", type="primary"):
            if guess:
                state.LockOut = True
                state.Lost = True
                state.FinalGuessm = guessm
                st.rerun()
            else:
                st.write("Oh come on you haven't even tried. Draw something and click Apply!")

        state.EndlessMW = st.slider("Endless Mode Molecular Weight Cutoff:", 100, 1000, int(state.EndlessMW))

# The Endless button
    with col2:
        if st.button(":green-background[♾️Endless Mode♾️]", type="secondary"):
            clean_slate()
            st.rerun()

# Similarity scoring
    if guess and not state.FirstEndless:
        tan = DataStructs.TanimotoSimilarity(targetfp, guessfp)
        HACdiff = targetHAC - guessHAC
        NumHDdiff = targetNumHD - guessNumHD
        NumHAdiff = targetNumHA - guessNumHA
        #RCdiff = targetRC - guessRC
        ArRCdiff = targetArRC - guessArRC
        AlRCdiff = targetAlRC - guessAlRC
 
        # check if guess is correct, output similarity stats if not
        if validguess and tan == 1:
            # TODO: splashier victory, write out copyable string
            #st.write("You got it! And it only took", state.guessnum, " tries!")
            #st.write("You got it! The answer is ",state.targetline[1].strip('\"'))
            state.LockOut = True
            state.Won = True
            state.FinalGuessm = guessm
            st.rerun()
        else:             
            st.write("Guess number:", state.guessnum)
            st.write("Your guess had this Tanimoto similarity to the target:", tan)
            st.write("Target minus guess heavy atom count:", HACdiff)
            st.write("Target minus guess H bond donors:", NumHDdiff)
            st.write("Target minus guess H bond acceptors:", NumHAdiff)
            #st.write("Target minus guess ring count:", RCdiff)
            st.write("Target minus guess aliphatic ring count:", AlRCdiff)
            st.write("Target minus guess aromatic ring count:", ArRCdiff)
            st.write("Atoms that share maximum common substructure with target highlighted below:")
 
        # If there was a validguess then update the table
        if validguess:
            # highlight max common substructure
            mcs_durl = view_mcs(targetm,guessm)
            st.image(mcs_durl)
            outrow = pd.DataFrame({"Guess Number": [int(state.guessnum)],
                                   "Tanimoto": [tan],
                                   "MCS": mcs_durl})
            state.outdf = pd.concat([state.outdf,outrow]) 
with col2:
    OutTable = st.dataframe(state.outdf,
                            column_config={"MCS": st.column_config.ImageColumn()},
                            hide_index=True)

# Just to deal with the guess autopopulation problem on Endless if the button is clicked 
state.FirstEndless = False
