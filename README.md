# SoLid visualisation code
This repository includes three codes used to visualize events for the SoLid antineutrino detector given the data incorperated in a specialised ROOT file. More information on the ROOT datasystem can be found on the website: https://root.cern/manual/storing_root_objects/. 

The first two programs, named resp. "writeToFile.py" and "visualFromFile.py" can be used to make a visualisation of a specific event, once the criteria for the event are met. The program is segmented in two parts in order to improve exportability; a user wishing to visualise a certain event, but prefers not to work in PYTHON, can still make use of the "visualFromFile.py"  code after creating a text file with the entry number of the event(s) which needed to be visualised. Furthermore in this repository is included the program "detectorVisualisation.ipynb" which can be used to make an animation of a period of detector time of the SoLid detector.

## writeToFile.py:
This program contains a function, named "WriteToFile", which writes a text file with all the events that meet the predefined requirements. This function takes a few arguments:

**Criteria:** This is a function able to check the criteria which an event need to meet for them to be visualised. This function should return a list with all the entries which meet the predefined criteria. The name of this function is left as a parameter so that, in order to visualise an event for which a function checking the criteria is not yet included in the file, a user only needs to include a function which returns a list of all events meeting the criteria in a given branch of the data to be able to use the WriteToFile function. A function checking for muon decay events resulting in a michel electron, named "michelCheck", and a fucntion checking for IBD events, named "IBDcheck", are already included in this program.

**filename:** This is the name that the eventual text file needs to have. This is left as a free parameter so that it can be optimized for the needed goals of a user without needing changes to the WriteToFile function.

**tree:** This is the name of the branch of the data tree that needs to be checked by the function "criteria". The data corresponding to different types of events is saved in different branches of the ROOT file. For example: an event which forms a cluster while passing through the detector is given in the branch "clusterTreeEntries" and an event where an electromagnetic signal (ES) and an nuclear signal (NS) are detected together, like an IBD event, is given in the branch "nsEsCoinsEntries". This parameter is already predefined for events from the "clusterTreeEntries" branch.

**numberOfEventsTogether:** This corresponds to the number of events that should be shown together. The code is implemented to be able to show either one or two events at the same time. In this way it is possible to show either a single particle, i.e. a muon, passing through the detector or a combination of two particles, i.e. an electron/positron and a neutron forming an IBD event or the decay of a muon resulting in a michel electron. 

## visualFromFile.py:
This program reads in the data within the text file produced by the "writeToFile.py" program and makes a visualisation of the event as well as the channels and the energy registered by each of the channels in projection. This function also takes in a few arguments:

**filename:** This is the name of the text file in which the data is saved, as produced by the "writeToFile.py" program. This is left as a free parameter so that a user that doesn't want to use the "writeToFile.py" program can still use the"visualFromFile.py" program as long as the data incorparated in the given text file is given in the same format as is produced by the "writeToFile.py" program.

**eventNumber:** This is the number of the event that needs to be shown. This way a user can choose which of the events within the text file should be visualised. 

**eventOneName:** This parameter corresponds to the name of the type of the first signal that should be shown. This is used to create a legend table in which all visualised events are incorperated.

**eventTwoName:** This parameter corresponds to the name of the type of a (potential) second signal that should be shown. This parameter is predefined so as to avoid errors when only one event should be shown.

**isMuon:** This is a boolian parameter that corresponds to whether or not a signal corresponds to a muon. If the event is indeed a muon, then the program will make a fit to the data. Therefore the needed fit parameters need to be requested from the data file

**NumberOfEventsTogether:**  (same as in the "writeToFile.py" program) This corresponds to the number of events that need to be shown together. This should correspond to the number of collumns in the text file.

## detectorVisualisation.ipynb:
The code included in the Jupyter Notebook file creates a interactive Jupyter Widget which can be used to visualize a time interval of data as collected by the SoLid detector. Each type of event is given a different color in order to distinguish them from each other; Blue for a muon, red for an electromagnetic signal (ES) and green for a nuclear signal (NS).
The widget creates two menu pages:

**Parameters:** Included in this menu page there are two parameter sliders. The first parameter, marked "time",  sets the duration of detector time that should be shown. Using this parameter it is possible to visualize all possible durations (in seconds) within the interval 0 ≲ t ≲ 5 s.The second parameter, marked "frames", makes it possible to change the amount of frames that should be used to make the animation. This is done to optimise the render time needed when making the animation. for example: In order to visualize a duration of up to  t = 1.8seconds there is a minimum of 1000 frames needed to make the animation. However When visualising 5 seconds of detector time, there are almost 3000 frames needed.

**Output:** In this menu page there are two dropdown menu's. The first dropdown menu, marked "Filled", allows the user to choose whether or not to just show the event as detected by the detector, or to fill the detector over time with all the events during the given detector time interval. The second dropdown menu, marked "Output", allows the user to choose whether or not to show the animation or to save it given the name given by the parameter "save name".

Once the parameters are chosen according to the correct specification of the user, the animation is made by pressing the "Create animation" button. 
