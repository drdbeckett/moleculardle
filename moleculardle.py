import streamlit as st
from streamlit_ketcher import st_ketcher
import streamlit.components.v1 as components
from streamlit_cookies_controller import CookieController
from datetime import datetime, timezone
from collections import defaultdict
from io import BytesIO
import pandas as pd
import pyperclip
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
from rdkit.Chem.Draw import rdMolDraw2D
rdDepictor.SetPreferCoordGen(True)

##############################################################
#     Moleculardle - A daily small molecule drawing game     #
#     Started 5/10/2024 - Developed by Daniel Beckett        #
##############################################################


state = st.session_state
st.set_page_config(layout='wide')
controller = CookieController()

# TODO: Change the colors automatically, maybe dark and purple?

# initialize some session state values ####################################
# just get the EST date
if 'today' not in state:
    local_tz = "US/Eastern"
    state.today = pytz.timezone(local_tz).localize(datetime.today())
    # Doing intial cookie parsing here to make sure we're not already playing or finished
    dateseed=str(state.today.day) + str(state.today.year) + str(state.today.month)
    GT_cookie = controller.get('GameToday')
    time.sleep(1)
    if GT_cookie:
    #    st.write("It stored: ", GT_cookie)
        # if we've not stored the cookie for today's game - reset cookies
        if GT_cookie != dateseed:
    #        st.write("Yesterday's news", GT_cookie)
            controller.set('GameToday',dateseed)
    else:
        controller.set('GameToday', dateseed)
    #    st.write("It didn't store")

    cookies = controller.getAll()
    #print("Final test: ", cookies)
    GT_cookie = controller.get('GameToday')
    #print("Final test: ", GT_cookie)

if 'hintguess' not in state:
    state.hintguess=""
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
    #TODO: Experiment with increasing readouts kept in table (ΔHBA, ΔHBD, ΔNAr, ΔNAl)
    state.outdf = pd.DataFrame({"Guess Number": [],
                                           "Tanimoto": [],
                                           "MCS": []})

st.title("MOLECULARDLE: A small molecule guessing game")
if not state.Endless:
    dailytitle="Daily Puzzle: "+str(state.today.month)+"/"+str(state.today.day)+"/"+str(state.today.year)
    st.title(dailytitle)
if state.Endless:
    st.title("🚨 ENDLESS MODE ACTIVATED 🚨")

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
# For visualizing the MCS and for getting the MCS hint
def view_mcs(targetm,guessm,hint):
    #rgba_color = (0.0, 0.0, 1.0, 0.2) # transparent blue
    colors = [(0.0, 0.0, 1.0, 0.2), (1.0, 0.0, 0.0, 0.2), (0.0, 1.0, 0.0, 0.3)] # transparent blue and red


    # We're doing 2 MCS calls, one loose then one strict
    loosemcs = rdFMCS.FindMCS([targetm,guessm], bondCompare=rdFMCS.BondCompare.CompareAny)
    loosemcs_mol = Chem.MolFromSmarts(loosemcs.smartsString)
    loose_gmatch = guessm.GetSubstructMatch(loosemcs_mol)
    tightmcs = rdFMCS.FindMCS([targetm,guessm], matchValences=True, ringMatchesRingOnly=True)
    tightmcs_mol = Chem.MolFromSmarts(tightmcs.smartsString)
    tight_gmatch = guessm.GetSubstructMatch(tightmcs_mol)
    tight_tmatch = targetm.GetSubstructMatch(tightmcs_mol)

    athighlights = defaultdict(list)
    arads = {}
    # These are for the hint
    keepatoms=[]
    #!keepbonds=[]

    # index of the match list is substructure index, entries are the atom indices
    # in either guessm or tightm, loop over the tight substructure match first
    for sid in range(len(tight_gmatch)):
        gid = tight_gmatch[sid]
        tid = tight_tmatch[sid]
        arads[gid] = 0.2
        # check if the guess aromaticity matches target
        #if guessm.GetAtomWithIdx(gid).IsInRing == True:
        #    if guessm.GetAtomWithIdx(gid).GetIsAromatic() == targetm.GetAtomWithIdx(tid).GetIsAromatic():
        #        athighlights[gid].append(colors[0])
        #        keepatoms.append(gid)
        #    else:
        #        athighlights[gid].append(colors[1])
        #else:
        #    athighlights[gid].append(colors[0])
        #    keepatoms.append(gid)
        if guessm.GetAtomWithIdx(gid).GetIsAromatic() == targetm.GetAtomWithIdx(tid).GetIsAromatic():
            athighlights[gid].append(colors[0])
            keepatoms.append(gid)
        else:
            athighlights[gid].append(colors[1])

    # Now check for anything from the loose substructure match
    for sid in range(len(loose_gmatch)):
        gid = loose_gmatch[sid]
        # check if the guess aromaticity matches target
        if gid not in athighlights:
            athighlights[gid].append(colors[1])
     
    #Now do the same procedure for bonds, loop over the tight bonds and get
    #indices we can use to get guess and target atoms involved in bonds and then
    #get those bond IDs
    bndhighlights = defaultdict(list)
    for bond in tightmcs_mol.GetBonds():
        gid1 = tight_gmatch[bond.GetBeginAtomIdx()]
        gid2 = tight_gmatch[bond.GetEndAtomIdx()]
        tid1 = tight_tmatch[bond.GetBeginAtomIdx()]
        tid2 = tight_tmatch[bond.GetEndAtomIdx()]
        gbid = (guessm.GetBondBetweenAtoms(gid1,gid2).GetIdx())
        tbid = (targetm.GetBondBetweenAtoms(tid1,tid2).GetIdx())
        # Check for aromaticity match
        #if guessm.GetBondWithIdx(gbid).IsInRing == True:
        #    if guessm.GetBondWithIdx(gbid).GetIsAromatic() == targetm.GetBondWithIdx(tbid).GetIsAromatic():
        #        bndhighlights[gbid].append(colors[0])
        #        #!keepbonds.append(gbid)
        #    else:
        #        bndhighlights[gbid].append(colors[1])
        #else:
        #    bndhighlights[gbid].append(colors[0])
        #    #!keepbonds.append(gbid)
        if guessm.GetBondWithIdx(gbid).GetIsAromatic() == targetm.GetBondWithIdx(tbid).GetIsAromatic():
            bndhighlights[gbid].append(colors[0])
            #!keepbonds.append(gbid)
        else:
            bndhighlights[gbid].append(colors[1])
     
    # Loose substructure match for bonds
    for bond in loosemcs_mol.GetBonds():
        gid1 = loose_gmatch[bond.GetBeginAtomIdx()]
        gid2 = loose_gmatch[bond.GetEndAtomIdx()]
        gbid = (guessm.GetBondBetweenAtoms(gid1,gid2).GetIdx())
        if gbid not in bndhighlights:
            bndhighlights[gbid].append(colors[1])
    
    #### Hint Gen! ####
    # Use the tight substructure search to generate the hint structure
    # Remove all atoms and bonds not highlighted blue
    # NOTE: originally had bond elimination in here but seems unnecessary, leaving for now with #! in front
    if hint:
        athighlights = defaultdict(list)
        bndhighlights = defaultdict(list)
        arads = {}
        badatoms=[]
        #!badbegin=[]
        #!badend=[]
        # Make list of atoms to remove
        for gatom in guessm.GetAtoms():
            gid = gatom.GetIdx()
            if gid not in keepatoms:
                badatoms.append(gid)        
        # Make lists of atoms in bonds to remove (between atoms that are kept)
        #!for gbond in guessm.GetBonds():
        #!    gbid = gbond.GetIdx()
        #!    if gbid not in keepbonds:
        #!        if gbond.GetBeginAtomIdx() in keepatoms:
        #!            if gbond.GetEndAtomIdx() in keepatoms:
        #!                badbegin.append(gbond.GetBeginAtomIdx())
        #!                badend.append(gbond.GetEndAtomIdx())

        badatoms.sort(reverse=True)

        # make editable version of guessm and remove bad atoms/bonds
        guessmw = Chem.RWMol(guessm)
        guessmw.BeginBatchEdit()
        #!for b1, b2 in zip(badbegin, badend):
        #!    guessmw.RemoveBond(b1,b2)
        for badatom in badatoms:
            guessmw.RemoveAtom(badatom)
        guessmw.CommitBatchEdit()
        if Chem.SanitizeMol(guessmw,catchErrors=True) == 0:
            guess = Chem.MolToSmiles(guessmw)
            guessm = Chem.MolFromSmiles(guess)
            state.guesses.append(guess)
            state.guessnum = state.guessnum+1
        else:
            return

        # Now we redo the MCS to get alignment with target again
        tightmcs = rdFMCS.FindMCS([targetm,guessm], matchValences=True, ringMatchesRingOnly=True)
        tightmcs_mol = Chem.MolFromSmarts(tightmcs.smartsString)
        tight_gmatch = guessm.GetSubstructMatch(tightmcs_mol)
        tight_tmatch = targetm.GetSubstructMatch(tightmcs_mol)
        # index of the match list is substructure index, entries are the atom indices
        # in either guessm or tightm, loop over the tight substructure match first
        for sid in range(len(tight_gmatch)):
            gid = tight_gmatch[sid]
            tid = tight_tmatch[sid]
            # if the number of bonds is less for an atom in guess then we highlight it
            if len(guessm.GetAtomWithIdx(gid).GetBonds()) < len(targetm.GetAtomWithIdx(tid).GetBonds()):
                athighlights[gid].append(colors[2])
                arads[gid] = 0.2
     
    # the actual draw command
    #mcs_pil = Draw.MolToImage(guessm, size=(400, 400), highlightAtoms=dict(athighlights), highlightBonds=dict(bndhighlights), highlightRadii=arads)
    d2d = rdMolDraw2D.MolDraw2DCairo(400,400)
    if hint:
        d2d.DrawMoleculeWithHighlights(guessm,"",dict(athighlights),dict(bndhighlights),arads,{})
   #     d2d.DrawMoleculeWithHighlights(guessm,"",{})
    else:
        d2d.DrawMoleculeWithHighlights(guessm,"",dict(athighlights),dict(bndhighlights),arads,{})
    #d2d.DrawMoleculeWithHighlights(guessm,"",dict(athighlights),dict(bndhighlights),arads,{})
    d2d.FinishDrawing()
    mcs_pil = BytesIO(d2d.GetDrawingText())

    # convert the png pil to a data url of a jpeg
    #buffered = BytesIO()
    #mcs_pil.save(buffered, format="JPEG")
    mcs_b64 = base64.b64encode(mcs_pil.getvalue()).decode("utf-8")
    #return 'data:image/jpeg;base64,' + mcs_b64
    return 'data:image/png;base64,' + mcs_b64

# For more human-oriented read out of guess propertie relative to target
def property_readout(num, prop):
    if num < -1:
        num = num*-1
        outstring="There are "+str(num)+" extra "+str(prop)+"s in the guess"
    elif num > 1:
        outstring="There are "+str(num)+" "+str(prop)+"s missing from the guess"
    elif num == 1:
        outstring="There is 1 "+str(prop)+" missing from the guess"
    elif num == -1:
        num = num*-1
        outstring="There is 1 extra "+str(prop)+" in the guess"
    elif num == 0:
        outstring="This guess has the correct number of "+str(prop)+"s!"
    return outstring

# For reading out specifically atom count difference (if different)
def atomcountdiffer(targetm,guessm,num,name):
    pstring = "[#"+str(num)+"]"
    namestring = name+" atom"
    pat = Chem.MolFromSmarts(pstring)
    delta = len(targetm.GetSubstructMatches(pat)) - len(guessm.GetSubstructMatches(pat))
    if delta != 0: st.write(property_readout(delta, namestring))
          
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
        elif int(i) == 9999:
            emojistring = emojistring+"🤫"
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
    state.guesses = True
    state.outdf = pd.DataFrame({"Guess Number": [],
                                "Tanimoto": [],
                                "MCS": []})

# For button color background
def ChangeButtonColour(widget_label, font_color, background_color='transparent'):
    htmlstr = f"""
        <script>
            var elements = window.parent.document.querySelectorAll('button');
            for (var i = 0; i < elements.length; ++i) {{ 
                if (elements[i].innerText == '{widget_label}') {{ 
                    elements[i].style.color ='{font_color}';
                    elements[i].style.background = '{background_color}'
                }}
            }}
        </script>
        """
    components.html(f"{htmlstr}", height=0, width=0)

## Exiting Function Town ##############################
#######################################################

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
# Canonicalize the smiles while reading it from target line    
target=Chem.CanonSmiles(state.targetline[3].strip('\"'))
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
              # 2. Today's working score and dataframe for the daily puzzle, JSON string?
              # 3. Win streak (maybe a lose streak and a daily streak later, keep it simple for the moment)
              # 4. Average number of guesses per win on Daily
              # 5. Average number of wins per Daily play
    if state.Won:
        st.write("You got it on guess ",str(state.guessnum),"! The answer is ",state.targetline[1].strip('\"'))
        link='https://go.drugbank.com/drugs/'+state.targetline[0].strip('\"')
        st.write(link)
        st.write(state.targetline[9].strip('\"'))
        st.image(view_mcs(state.FinalGuessm,targetm,False))
        if not state.Endless:  
            st.write("Copy the emoji string to show off to your friends, colleagues, and enemies!")
            emojistring = emojify()
            st.write("Moleculardle ",str(state.today.month),"/",str(state.today.day),"/",str(state.today.year),": ", emojistring)
            outstring = "Moleculardle "+str(state.today.month)+"/"+str(state.today.day)+"/"+str(state.today.year)+": "\
                    +emojistring+"\n I got it in "+str(state.guessnum)+\
                    "! Think you can do better? http://moleculardle.streamlit.app"
            if st.button("Copy to Clipboard", type="secondary"):
                pyperclip.copy(outstring)

    if state.Lost:
        # TODO: Move the image up and try to make it larger
        st.write("Better luck next time! Here's how close you got.")
        st.write("The answer is ",state.targetline[1].strip('\"'))
        link='https://go.drugbank.com/drugs/'+state.targetline[0].strip('\"')
        st.write(link)
        st.write(state.targetline[9].strip('\"'))
        st.image(view_mcs(state.FinalGuessm,targetm,False))
        if not state.Endless:  
            st.write("Copy the emoji string to demonstrate how hard you tried before tapping out!")
            emojistring = emojify()
            st.write("Moleculardle ",str(state.today.month),"/",str(state.today.day),"/",str(state.today.year),": ", emojistring)
            outstring = "Moleculardle "+str(state.today.month)+"/"+str(state.today.day)+"/"+str(state.today.year)+": "\
                    +emojistring+"\n I gave in after "+str(state.guessnum)+\
                    "! Think you can do better? http://moleculardle.streamlit.app"
            if st.button("Copy to Clipboard", type="secondary"):
                pyperclip.copy(outstring)

    # reset button that initiates endless mode
    # TODO: Make this and the other endless button primary buttons that are green
    if st.button("🚨 Continue in Endless Mode? 🚨", type="secondary"):
        clean_slate()
        st.rerun()

    ChangeButtonColour('🚨 Continue in Endless Mode? 🚨', 'white', 'green')
    state.EndlessMW = st.slider("Endless Mode Molecular Weight Cutoff", 100, 1000, int(state.EndlessMW))
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
        if st.button("☠️  Give up? ☠️", type="primary"):
            if guess:
                state.LockOut = True
                state.Lost = True
                state.FinalGuessm = guessm
                st.rerun()
            else:
                st.write("Oh come on you haven't even tried. Draw something and click Apply!")
        
        # The Hint button
        if guess:
            if st.button("🤫  Need a hint? 🤫", type="secondary"):
                hint = view_mcs(targetm,guessm,True)
                if hint:
                    tan = 9999
                    outrow = pd.DataFrame({"Guess Number": [int(state.guessnum)],
                                           "Tanimoto": [tan],
                                           "MCS": hint})
                    state.outdf = pd.concat([state.outdf,outrow])
                    st.markdown('''This is the strictest maximum common substructure of your last guess
                                   with atoms that need to be grown from highlighted with green circles.
                                   This image will populate into the table and costs you 1 guess. Good luck!''')
                    st.image(hint)
                else:
                    st.markdown('''Ran into an issue decomposing the molecule! No hint for you this time -
                                   odds are your guess has an aromatic ring that is broken when the blue
                                   substructure is considered. Fix this and try again, this did not cost you a guess.''')

# The Endless molecular weight slider and button
        state.EndlessMW = st.slider("Endless Mode Molecular Weight Cutoff:", 100, 1000, int(state.EndlessMW))

        if st.button("🚨 Initiate Endless Mode? 🚨", type="secondary"):
            clean_slate()
            st.rerun()

    ChangeButtonColour('🚨 Initiate Endless Mode? 🚨', 'white', 'green')

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
            # TODO: splashier victory
            #st.write("You got it! And it only took", state.guessnum, " tries!")
            #st.write("You got it! The answer is ",state.targetline[1].strip('\"'))
            state.LockOut = True
            state.Won = True
            state.FinalGuessm = guessm
            st.rerun()
        else:             
            st.write("Guess number:", str(state.guessnum))
            with st.popover("Empirical formula comparison"):
                guessformula = rdMolDescriptors.CalcMolFormula(guessm)
                st.write("Guess empirical formula: ",guessformula)
                if guessformula == targetformula:
                    st.write("You've nailed the formula! Now you've just got to get the rest right.")
                else:  
                    atomcountdiffer(targetm,guessm,6,"carbon")
                    atomcountdiffer(Chem.AddHs(targetm),Chem.AddHs(guessm),1,"hydrogen")
                    atomcountdiffer(targetm,guessm,7,"nitrogen")
                    atomcountdiffer(targetm,guessm,8,"oxygen")
                    atomcountdiffer(targetm,guessm,9,"fluorine")
                    atomcountdiffer(targetm,guessm,15,"phosphorus")
                    atomcountdiffer(targetm,guessm,16,"sulfur")
                    atomcountdiffer(targetm,guessm,17,"chlorine")
                    atomcountdiffer(targetm,guessm,35,"bromine")
              
            st.write("Your guess had this Tanimoto similarity to the target:", str(tan))
            st.write(property_readout(HACdiff, "heavy atom"))
            st.write(property_readout(NumHDdiff, "H-bond donor"))
            st.write(property_readout(NumHAdiff, "H-bond acceptor"))
            #st.write(property_readout(RCdiff, "ring"))
            st.write(property_readout(AlRCdiff, "*aliphatic* ring"))
            st.write(property_readout(ArRCdiff, "*aromatic* ring"))
            substruct_line = '''Atoms/bonds in the maximum common substructure with the target molecule are highlighted below,
                       atoms/bonds highlighted in  :blue[blue] matched a :blue[tight] substructure search while atoms/bonds
                       in :red[red] only matched a :red[loose] substructure search and need some work.'''
            st.markdown(substruct_line)
 
        # If there was a validguess then update the table
        if validguess:
            # highlight max common substructure
            mcs_durl = view_mcs(targetm,guessm,False)
            st.image(mcs_durl)
            outrow = pd.DataFrame({"Guess Number": [int(state.guessnum)],
                                   "Tanimoto": [tan],
                                   "MCS": mcs_durl})
            state.outdf = pd.concat([state.outdf,outrow])
              
        # Substructure highlighting help guide
        with st.popover("Substructure highlighting guide"):
            st.markdown(
                    """
                    There are several key differences between the :blue[blue] and :red[red] highlighting,
                    in all cases :blue[blue] is a stricter substructure search and overrides :red[red]. Below
                    are some basic guidelines with example highlightings overlaid on guess and target molecules.
                    - The valence of atoms can be wrong in :red[red], as well as the bond order, but not in :blue[blue].
                    """
                    )
            st.image("examples/ex1.png")
            st.markdown(
                    """
                    - If the target is a ring but guess is not a ring, :red[red] will match but :blue[blue] will not.
                    """
                    )
            st.image("examples/ex2.png")
            st.markdown(
                    """
                    - If the target ring is aliphatic but the guess is aromatic, :red[red] will match but :blue[blue] will not (and vice versa).
                    """
                    )
            st.image("examples/ex3.png")
            st.markdown(
                    """
                    - If the guess ring is not the right size then :blue[blue] will match until it reaches the closing bond (for smaller rings)
                    or until it reaches the final atom of the target ring (for larger rings). :red[Red] will try to find the longest bond path
                    around the full molecule and may highlight bonds :blue[blue] does not. This can lead to behavior that is not obvious at first
                    glance, especially in complicated cases where a ring fusion is in the incorrect position. If the daily challenge is too difficult
                    it can be helpful to play on Endless with a lower molecular weight to hone your skills. I believe in you!
                    """
                    )
            st.image("examples/ex4.png")
              
with col2:
    OutTable = st.dataframe(state.outdf,
                            column_config={"MCS": st.column_config.ImageColumn()},
                            hide_index=True)

state.FirstEndless = False
