import streamlit as st
from streamlit_ketcher import st_ketcher
from datetime import date
from io import BytesIO
import pandas as pd
import linecache
import random
import base64
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
# initialize some session state values ####################################
if 'guessnum' not in state:
    state.guessnum = 0

# Block for Win/Loss condition    
if 'LockOut' not in state:
    state.LockOut = False
if 'Won' not in state:
    state.Won = False
if 'Lost' not in state:
    state.Lost = False

# store list of guessed SMILES strings to quickly check if a guess is unique
if 'guesses' not in state:
    state.guesses = []
if 'outdf' not in state:
    state.outdf = pd.DataFrame({"Guess Number": [],
                                           "Tanimoto": [],
                                           "MCS": []})
# Generate the line number for grabbing the target from the current date     
# and pull out the full line as an array to parse at our leisure
# 0 = DBID, 1 = name, 3 = SMILES, 9 = summary
if 'targetline' not in state:
    # There's probably a better way than hardcoding number of lines but also probably fastest
    numlines=8288
    dateseed=str(date.today().day) + str(date.today().year) + str(date.today().month)
    random.seed(dateseed)
    targetnum=random.randrange(1,numlines)
    #tline = linecache.getline('/Users/dbeckett/tutorials/dsmiles_cleaned.csv', targetnum) 
    tline = linecache.getline('dsmiles_cleaned.csv', targetnum) 
    state.targetline = [ '"{}"'.format(x) for x in list(csv.reader([tline], delimiter=',', quotechar='"'))[0] ]
    # debug block 
    #print(targetline[1])
    #print(targetline[3])
    #print(targetline[9])
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

    # convert the png pil to a data url of a jpeg
    buffered = BytesIO()
    mcs_pil.save(buffered, format="JPEG")
    mcs_b64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return 'data:image/jpeg;base64,' + mcs_b64

def emojify():
    emojistring=""
    for i in state.outdf['Tanimoto']:
     #   st.write(i['Tanimoto'])
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


#######################################################
# Target definition and initialization
### debug targets
#target = "NC1CCCCC1"
#target = "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4"
#target = "O=C(N(C)C1=C2N(C)C=N1)N(C)C2=O"

# Get the target SMILES string from the csv file line
# usually have to strip the double quotes
target=state.targetline[3].strip('\"')
targetm = Chem.MolFromSmiles(target)
targetHAC = rdMolDescriptors.CalcNumHeavyAtoms(targetm)
targetNumHD = rdMolDescriptors.CalcNumHBD(targetm)
targetNumHA = rdMolDescriptors.CalcNumHBA(targetm)
targetRC = rdMolDescriptors.CalcNumRings(targetm)
targetfp = FingerprintMols.FingerprintMol(targetm)
targetformula = rdMolDescriptors.CalcMolFormula(targetm)


guess = ""
validguess = False

# for the winners or losers
if state.Won:
    st.write("You got it on guess ",str(state.guessnum),"! The answer is ",state.targetline[1].strip('\"'))
    link='https://go.drugbank.com/drugs/'+state.targetline[0].strip('\"')
    st.write(link)
    st.write(state.targetline[9].strip('\"'))
    st.image(view_mcs(state.FinalGuessm,targetm))
    st.write("Copy the emoji string to show off to your friends, colleagues, and enemies!")
    st.write("Structurdle ",str(date.today().month),"/",str(date.today().day),"/",str(date.today().year),": ", emojify())

if state.Lost:
    st.write("Better luck next time! Here's how close you got.")
    st.write("The answer is ",state.targetline[1].strip('\"'))
    link='https://go.drugbank.com/drugs/'+state.targetline[0].strip('\"')
    st.write(link)
    st.write(state.targetline[9].strip('\"'))
    st.image(view_mcs(state.FinalGuessm,targetm))
    st.write("Copy the emoji string to demonstrate how hard you tried before tapping out!")
    st.write("Structurdle ",str(date.today().month),"/",str(date.today().day),"/",str(date.today().year),": ", emojify())
###########################

# Get input structure and properties
# TODO: style these writes
if not state.LockOut:
    st.write("Guess the drug!")
    st.write("Target empirical formula:", targetformula)
    guess = st_ketcher()

# get properties from the input guess
if guess: 
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
    guessRC = rdMolDescriptors.CalcNumRings(guessm)

# The give up button
if not state.LockOut:
    if st.button("Give up?", type="primary"):
        if guess:
            state.LockOut = True
            state.Lost = True
            state.FinalGuessm = guessm
            st.rerun()
        else:
            st.write("Oh come on you haven't even tried. Draw something and click Apply!")

# Similarity scoring
    if guess:
        tan = DataStructs.TanimotoSimilarity(targetfp, guessfp)
        HACdiff = targetHAC - guessHAC
        NumHDdiff = targetNumHD - guessNumHD
        NumHAdiff = targetNumHA - guessNumHA
        RCdiff = targetRC - guessRC
 
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
            st.write("Target minus guess ring count:", RCdiff)
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
OutTable = st.dataframe(state.outdf,
                        column_config={"MCS": st.column_config.ImageColumn()},
                        hide_index=True)

