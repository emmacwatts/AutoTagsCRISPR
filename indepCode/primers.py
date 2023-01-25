
# CONDITIONS VARIABLES - [stringent[min,opt,max], relaxed, desperate]
#primer_position = ["start(ATG)","stop"]

primer_name = ["val-F","HAL-F","HAL-R","HAR-F","HAR-R","val-R"]
initial_primer_region = [[0,100],[200,200],[1270,29],[1303,30],[2200,200],[2500,100]]
extended_primer_region = [[0,100],[200,300],[1270,29],[1303,30],[2200,200],[2500,100]]

GC_content = [[30,50,70],[25,50,75],[20,50,80]]
size = [[18,20,25],[18,20,26],[17,20,30]]
max_end_GC = [3,4,4]
GC_clamp = [1,1,0]
TH_max_hairpin = [47.00,48.00,72.00]
max_polyx = [5,6,8] 
stringency_levels = 3

def DesignPrimer(template,stringency,GC_content,max_end_GC,size,GC_clamp,TH_max_hairpin,max_polyx,
                 primer_region, primer_number):
 
  '''
  THIS FUNCTION CALCULATES A PRIMER GIVEN THE CONDITIONS AND HOW STRICT
   THOSE CONDITIONS SHOULD BE

  stringency is a list with n elements (i.e., n different stringency levels). Each condition is a 
  list (or a list of list). The list contains n elements (one for each striengency level). The 
  innermost list in a list of list contains [min,opt,max]
  
  '''
  #I've redefined primer type this way throughout to cut some code
  primer_type = primer_name[-1]
  if primer_type == "F":
    left = 1
    right = 0
  elif primer_type == "R":
    left = 0
    right = 1

  primer = p3.bindings.designPrimers(
    {
        'SEQUENCE_TEMPLATE': template,
        'SEQUENCE_INCLUDED_REGION': primer_region, 
    },
    {
        'PRIMER_NUM_RETURN':primer_number,

        'PRIMER_TASK': "generic",
        'PRIMER_PICK_LEFT_PRIMER': left,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_PICK_RIGHT_PRIMER': right,

        'PRIMER_MIN_GC': GC_content[stringency][0],
        'PRIMER_OPT_GC_PERCENT': GC_content[stringency][1],
        'PRIMER_MAX_GC': GC_content[stringency][2],
     
        'PRIMER_MIN_SIZE': size[stringency][0],
        'PRIMER_OPT_SIZE': size[stringency][1],
        'PRIMER_MAX_SIZE': size[stringency][2],
     
        'PRIMER_MAX_END_GC': max_end_GC[stringency],
     
        'PRIMER_GC_CLAMP': GC_clamp[stringency],    
             
        'PRIMER_MAX_HAIRPIN_TH':TH_max_hairpin[stringency],
    
        'PRIMER_MAX_POLY_X': max_polyx[stringency],
    })
  
  return primer

def primers_cleanup_dict(Gene_ID = 'Gene_ID', position = 'start_stop', primer_type = 'primer_name', primer_sequence = 'primer_sequence', stringency_level = 'stringency'):
    primers_cleanup = {'Gene_ID': Gene_ID, 'position': position, 'primer_type': primer_type, 'primer_sequence': primer_sequence, 'stringency_level': stringency_level}
    return primers_cleanup

def DoIHaveAPrimer(extension,Gene_ID,start_stop,primer_type, primer_name,stringency,primer):
  
  if extension == 0:
    Extended_tag = ""
  elif extension == 1:
    Extended_tag = "e"

  if primer_type == "F":
    if primer['PRIMER_LEFT_NUM_RETURNED']>0:
      primers_cleanup = primers_cleanup_dict(primer_sequence = primer['PRIMER_LEFT_0_SEQUENCE'], stringency_level = f'{Extended_tag}{stringency+1}')
      warning_variable = False
    else:
      primers_cleanup = ()
      warning_variable = True

  elif primer_type == "R":
    if primer['PRIMER_RIGHT_NUM_RETURNED']>0:
        primers_cleanup = primers_cleanup_dict(primer_sequence = primer['PRIMER_RIGHT_0_SEQUENCE'], stringency_level = f'{Extended_tag}{stringency+1}')
        warning_variable = False
    else:
      primers_cleanup = ()
      warning_variable = True

  return primers_cleanup, warning_variable 

def SixPrimersCalculator(Gene_ID,start_stop,template,primer_name,initial_primer_region,enlarged_primer_region,GC_content,size,max_end_GC,GC_clamp,TH_max_hairpin,max_polyx,stringency_levels):
  primers_table = pd.DataFrame() # DataFrame used to save the 6 primers
  
  for pt in range(0,len(primer_name)): #loop to get as many primers as specified in the primer_name list (here 6)
    
    primer_type = primer_name[-1]

    # defining how many primer outputs I want to have (one is NOT enough for the special cases)
    if ((pt==2) or (pt==3)):  # output is only one primer if we do not care about the exact initial position (i.e., for 4/6 primers)
      primer_number = 1000
    else:
      primer_number = 1

   ###############################################################################
    # loop to ACTUALLY take care of the primers (one at a time, considering different stringencies)
    stringency = 0 #initial condition (most stringent conditions)

    safety_net_primer_forward = ()  # empty variables to save HAL-R and HAR-F primers at a random position (BUT still in the right range)
    safety_net_primer_reverse = ()  # this is useful in case we can't find any satisfactory primers in the right position
    primer_region = initial_primer_region

    while stringency<stringency_levels: # TO DO NEXT: must add extra level to go in larger dimensions (with if statement to change the borders soon below)

      primer = DesignPrimer(template,stringency,GC_content,max_end_GC,size,GC_clamp,TH_max_hairpin,max_polyx,    # designing the primer
                          primer_region[pt], primer_type, primer_number)

      if primer_type == "F": # necessary only for special cases 2 and 3 - the # of calculated primers allows us to loop over them and choose the best
        number_of_calculated_primers = primer['PRIMER_LEFT_NUM_RETURNED']  
      elif primer_type == "R":
        number_of_calculated_primers = primer['PRIMER_RIGHT_NUM_RETURNED']
          
      ####### special cases -  HAL-R and HAR -F #############
      if (pt == 2 or pt == 3): # TO DO NEXT: when we add the extra level this will need to be skipped --> must add and stirngency<stringency level
        primers_cleanup = () # empty variable to save the final primer of choice

        for sc in range(0,number_of_calculated_primers): # evaluating all the primers obtained

            if (pt==2): # HAL-R

              # saving a random "good" primer as a "safety net" in case we can not find any other in the exact position
              if safety_net_primer_reverse == ():   # saving of a "net" happens only if we find a primer + no primers can get overwritten
                safety_net_primer_reverse = primer['PRIMER_RIGHT_0_SEQUENCE'] # saving the first one - for each primer the safety net is overwritten
              
              temporary_name_string = f'PRIMER_RIGHT_{sc}_SEQUENCE'  
              reverse_complement = revComp(primer[temporary_name_string]) # checking position of the primer
              primer_position = template.index(reverse_complement)
              if (primer_position + len(reverse_complement)) == 1300: 
                primers_cleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = primer[f'PRIMER_RIGHT_{sc}_SEQUENCE'], stringency_level = stringency+1)

            if (pt==3): #HAR-F

              # saving a random "good" primer as a "safety net" in case we can not find any other in the exact position
              if safety_net_primer_forward == ():
                safety_net_primer_forward = primer['PRIMER_LEFT_0_SEQUENCE']
            
              temporary_name_string = f'PRIMER_LEFT_{sc}_SEQUENCE'
              primer_position = template.index(primer[temporary_name_string])
              if primer_position == 1303: # checking primer position
                primers_cleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = primer[f'PRIMER_LEFT_{sc}_SEQUENCE'], stringency_level = stringency+1)

        # TO DO: careful with stringency levels after extenstion
        if stringency == (stringency_levels-1) and primers_cleanup == (): # i.e. if we are at the last stringency and we still haven't left the loop....
          if (pt==2):
            if safety_net_primer_reverse == ():
              last_resort = template[1275:1300]
              primer_last_resort = revComp(last_resort)
              primers_cleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = primer_last_resort, stringency_level = "25bp no conditions")

            else:
              safety_net_primer_reverse_t = revComp(safety_net_primer_reverse)
              net_position = template.index(safety_net_primer_reverse_t)
              extened_safety_net_rev_comp = template[net_position:1300]
              extended_safety_net_primer_reverse = revComp(extened_safety_net_rev_comp)
              primers_cleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = extended_safety_net_primer_reverse, stringency_level = f"extended from {stringency+1}")

          if (pt==3):
            if safety_net_primer_forward == ():
                primers_cleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = template[1303:1328], stringency_level = "25bp no conditions")
            else:
              net_position = template.index(safety_net_primer_forward)
              extension = template [1303:net_position]
              extended_safety_net_primer_forward = extension+safety_net_primer_forward
              primers_cleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = extended_safety_net_primer_forward, stringency_level = f"extendedb from {stringency+1}")

        if primers_cleanup == ():
          stringency +=1
        else:
          primers_cleanup_table = pd.DataFrame([primers_cleanup]) # coverting to dataframe
          primers_table = pd.concat([primers_table,primers_cleanup_table]) # extending existing list with dataframe that was just created
          stringency = stringency_levels # TO DO: SEE HOW TO MODIFY THIS NEXT

      ########### normal cases #########
      if (pt == 0 or pt == 1 or pt == 4 or pt == 5):
        if primer_region == initial_primer_region:
          extension = 0
        else:
          extension = 1

        primers_cleanup, warning_variable = DoIHaveAPrimer(extension,Gene_ID,start_stop,primer_type, primer_name[pt],stringency,primer)      # function 2 created by me 
        
        if warning_variable == False: # i.e., if I do not get a warning then we are done --> we leave the while loop with 
          primers_cleanup_table = pd.DataFrame([primers_cleanup]) # coverting to dataframe
          primers_table = pd.concat([primers_table,primers_cleanup_table]) # extending existing list with dataframe that was just created
          stringency = stringency_levels # leaving the loop
        else:
          if stringency == (stringency_levels-1):
            if primer_region == enlarged_primer_region:
                primerscleanup = primers_cleanup_dict(primer_type = primer_name[pt], primer_sequence = "primer could not be calculated",stringency_level = "NA")
                stringency = stringency_levels # leaving the loop

            elif primer_region == initial_primer_region:
              stringency = 0  
              primer_region = enlarged_primer_region
          else:
            stringency +=1

  return(primers_table)    
