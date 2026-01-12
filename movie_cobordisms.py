import copy


# HELPER FUNCTIONS

def canonical_perm(n, mappings):

    """
    Given a list of mappings, recovers a permutation {x:y} such that the mappings are satisfied by the permutation. 

    INPUTS:
    n: Input length (permutations of {0,1,ldots,n-1}
    mappings: A dictionnary of the form {x:y} where we desire that x is mapped to y

    OUTPUTS:
    A dictionnary {x:y} representing a permutation where x is mapped to y

    REMARK: This assumes the pairing in mappings is injective.
    """
    permutation = {} 
    keys,values = mappings.keys(), mappings.values()
    l0,l1 = [],[] 
    
    for i in range(n):
        if i not in keys:
            l0.append(i)
        if i not in values:
            l1.append(i)
            
    if len(l0) != len(l1):
        raise ValueError("mapping given is not injective")

    for index in range(len(l0)):
        inp,out = l0[index],l1[index]
        if inp!= out:
            permutation[inp] = out 

    for k in keys:
        permutation[k]= mappings[k] 

    return permutation
    
    

def strand_index(crossing, strand, direction):
    
    """
    Given a crossing in PD notation and a strand in that crossing, it will
    return the index (0,1,2 or 3) of the strand in the crossing. 
    This assumes the strand is in the crossing.
    If direction == 1, it will return the index at which the strand is incoming.
    Else, it will return the index at which the strand is outgoing.

    INPUTS:
    strand: The strand we want to check the direction of
    crossing: The desired crossing as a list of 4 strands 
    direction: 1 for incoming, -1 for outgoing

    OUTPUTS:
    The index of the strand at the crossing for which it is incoming/outgoing.
    
    Eg:
    strand_index([1,1,2,2], 1, 1) = 0, strand_index([1,1,2,2],1,-1) = 1
    """
    
    matches = len([x for x in crossing if x == strand])
    if matches == 1:
        return crossing.index(strand)

    if direction == 1:
        # incoming strand 
        if crossing[0] == strand:
            return 0
        return (1 if crossing[1] == strand else 3)
    else: 
        # outgoing strand
        if crossing[2] == strand:
            return 2 
        return (1 if crossing[1] == strand else 3)

def strand_circle(strand, state):

    """
    Given an enhanced Kauffman state, return the sign of the circle having "strand" in its boundary.
    WARNING: This assumes the strand indeed exists in the state's notation. Furthermore, if state is a state of the unknot, it will return
    the sign of the unknot's unique circle.

    INPUTS:
    strand: The strand number in PD notation.
    state: The enhanced Kauffman state.

    OUTPUTS:
    A tuple (c,sign) where c is the circle containing the strand in the enhanced Kauffman state, 
    and sign is the sign of that circle.
    
    """

    if not len(state[0]):
        return (tuple(),-1) if len(state[1]) else (tuple(),1)

    for circle in state[1]:
        for branch in circle:
            if strand in branch:
                return (circle,-1)

    for circle in state[2]:
        for branch in circle:
            if strand in branch:
                return (circle,1)

    raise ValueError("Branch not in state")


def strand_sign(strand, state):

    """
    Given an enhanced Kauffman state, return the sign of the circle having "strand" in its boundary.
    WARNING: This assumes the strand indeed exists in the state's notation. Furthermore, if state is a state of the unknot, it will return
    the sign of the unknot's unique circle.

    INPUTS:
    strand: The strand number in PD notation.
    state: The enhanced Kauffman state.

    OUTPUTS:
    1 for a positive circle containing the strand, -1 for a negative one.  
    """

    return strand_circle(strand, state)[1]

def delete_circles(state, strands):

    """
    Given an enhanced Kauffman state and list of strands, returns the circles not containing said strands.

    INPUTS:
    state: an enhanced Kauffman state
    strands: A list of strands

    OUTPUTS:
    A list of two sets [negative_circles, positive_circles] where each circle does not contain the mentioned strands.
    """
    negative_circles = set()
    positive_circles = set()
    strands = set(strands)
    for neg_circle in state[1]:
        indicator = True
        for branch in neg_circle:
            if set(branch).intersection(strands):
                indicator = False
        if indicator:
            negative_circles.add(neg_circle)

    for pos_circle in state[2]:
        indicator = True
        for branch in pos_circle:
            if set(branch).intersection(strands):
                indicator = False
        if indicator:
            positive_circles.add(pos_circle)  

    return [negative_circles, positive_circles]
        



# MOVIE CLASS

class Movie():
    
    """
    Class representing movie cobordisms. Explanation of attributes below:

    
    -links : The links of the movie in order as a list.
    -current_min_index : The minimum index used in the PD notations for the last link
    -plot_index : Internal tracker for the index of a plot (refer to next_link() method)
    -last_degree : Movie moves are graded by Euler characteristic. Starting with a certain quantum degree, this will change throughout movie moves. This gives the last qdegree.
    -maps : Chain maps (maps[i] is the chainmap from the i'th complex to the i+1'st complex)
    -ring : Underlying ring for Khovanov homology
    -bases : Bases for the Khovanov homology (bases[i] gives a dictionnary who'se keys are (homological_deg, quantum_deg) and who'se values are the ordered enhanced Kauffman states)
    -complexes : List of chain complexes of the links at the appropriate qdegree.

    Explanation of methods:
    -twist : Does an R1 move on a strand.
    -untwist : Undo an R1 move 
    -poke : Does an R2 move
    -unpoke : Undo an R2 move 
    -slide : Does an R3 move 
    -birth : Birthes a twist 
    -death : Kills a twist 
    -saddle : Does a saddle move on a strand
    -next_link : Outputs the (absolute) next link. By absolute, I mean that all PD notation's signs are taken in absolute values for plotting purposes.
    -reset_plot_index : Sets plot_index back to 0.
    
    See methods documentation for further reference.

    Example of usage: (Trefoil movie)
    
    L = Link([])
    
    my_movie = Movie(L, starting_qdeg=1)
    my_movie.twist(-1,-1,-1).pd_code()
    
    # Output: [[-1, -2, -2, -1]]
    
    my_movie.twist(-1,-1,-1).pd_code()
    
    # Output:  [[-4, -2, -2, -1], [-1, -3, -3, -4]]
    
    my_movie.twist(-4,-1,-1).pd_code()
    
    # Output: [[-6, -2, -2, -1], [-1, -3, -3, -4], [-4, -5, -5, -6]]
    
    my_movie.saddle(-2,-3).pd_code()
    
    # Output: [[-6, -3, -2, -1], [-1, -2, -3, -4], [-4, -5, -5, -6]]
    
    my_movie.saddle(-3,-5).pd_code()
    
    # Output: [[-6, -5, -2, -1], [-1, -2, -3, -4], [-4, -3, -5, -6]]

    

    Repeatedly running movie.next_link().plot() will plot the movie (in Jupyter for example)

    Remark: We take the convention of using negative integers for strands in PD notation.
    """

    # DEFINITION

    def __init__(self, link, starting_qdeg = 1, ring=QQ):
        
        """
        Initializes the movie.
        
        INPUTS : 
        link : Starting link for your movie.
        starting_qdeg : Starting quantum degree for which we compute the chain maps on khovanov homology (1 by default)
        ring : Underlying ring for Khovanov homology
        """
        
        self.links = [link]  #
        self.current_min_index = 0 if not link.pd_code() else min([x for crossing in link.pd_code() for x in crossing]) #
        self.plot_index = 0
        self.last_degree = starting_qdeg #
        self.maps = [] #
        self.ring = ring 
        self.bases = [] #
        self.complexes = [] #

        initial_complex, initial_base = height_khovanov_chain_complex(link, starting_qdeg, True, ring=self.ring)
        self.bases.append(initial_base)
        self.complexes.append(initial_complex)

    # REIDEMEISTER MOVES

    # R1 
    
    def twist(self, strand, orientation = 1, strand_type = 1, print_pd=True):
        
        """
        Creates a twist (Reidemeister 1 Move)
        Computes the new induced link as well as the appropriate chain maps.
        
        INPUTS: 
        strand: Strand number of the last link in PD notation
        orientation : If we want the twist's loop's 0-resolution to contain a loop or not. If orientation = 1, it will, otherwist it won't. 1 by default.
        strand_type : If we want the loop's crossing to be overstrand (strand_type = 1) or understrand (strand_type = -1)
        print_pd : If we want to output the pd code of the final link or not
    
        OUTPUTS: 
        Returns the last link.


        REMARK : The convention for number the new strand is as follows: If the strand we want to loop goes from crossing i to crossing j, then
        after the loop, the strand outgoing from i remains with the same label. The loop becomes the minimum label-1 and the strand going to crossing j
        becomes with minimum label -2. The minimum label is then updated.
        """

        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)
        last_degree = self.last_degree
        
        twist_knot_crossing = []
        loop_label = -2 if not crossings else self.current_min_index - 1
        strand = -1 if not crossings else strand
        new_strand_label = -1 if not crossings else loop_label - 1
        self.current_min_index = -2 if not crossings else new_strand_label

        if orientation == 1 and strand_type == 1 :
            twist_knot_crossing = [loop_label, loop_label, new_strand_label, strand]
        elif orientation == 1 and strand_type == -1:
            twist_knot_crossing = [strand, new_strand_label, loop_label, loop_label]
        elif orientation == -1 and strand_type == 1:
            twist_knot_crossing = [loop_label, strand, new_strand_label, loop_label]
        else :
            twist_knot_crossing = [strand, loop_label, loop_label, new_strand_label]

        new_crossings.append(twist_knot_crossing)

        # Computing new pd notation

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else len(orientationnes) - n_minus
        
        if crossings:
        
            directions = last_link._directions_of_edges()
            if strand not in directions[0]:
                raise ValueError("Invalid input provided: The strand is not in the last link's PD notation.")
                
            incoming_to = directions[1][strand]
            j = crossings.index(incoming_to)
            new_crossings[j][strand_index(new_crossings[j], strand, 1)] = new_strand_label

        new_link = Link(new_crossings)
        self.links.append(new_link)

        # Computing chain map
        
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1] # arrange them by (i,j)
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree, True, ring=self.ring) # arrange them by (i,j)
        self.complexes.append(final_complex)
        self.bases.append(final_bases)

        chain_maps = {}
         
        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree) not in final_bases) else len(final_bases[(homological_degree, last_degree)]) 
            m = matrix(self.ring, domain_size, image_size)
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree)][image_state_index]

                    if image_state[0][:-1] == state[0]:

                        remaining_circles1 = delete_circles(state, [strand])
                        remaining_circles2 = delete_circles(image_state, [strand, loop_label, new_strand_label])
                        
                        if orientation == 1:
                            if image_state[0][-1] == 0 and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])):
                                m[state_index,image_state_index] = -1 if (((loop_label, loop_label),) in image_state[2]) else 1
                        else:
                            if image_state[0][-1] == 1 and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and (((loop_label, loop_label),) in image_state[2]):
                                m[state_index,image_state_index] = 1 

            chain_maps[homological_degree]= m.transpose()

        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)

        if print_pd:
            print(new_crossings)
        return new_link

    def untwist(self, loop_label, print_pd = True):
        
        """
        Untwists a link, opposite to the twist (R1) move. 

        INPUTS:
        loop_label: The label of the loop strand.
        print_pd : If we want to output the pd code of the final link or not
        
        OUTPUTS:
        Returns the last link.

        WARNING: If you untwist a twisted unknot that is not the entire link, an error will be thrown. Refer to the death method instead in this 
        specific instance.
        """
        
        # Checking input validity

        bad_crossings = self.links[-1].pd_code()
        ncrossings = len(bad_crossings)

        directions = self.links[-1]._directions_of_edges()
        loop_crossing = directions[0][loop_label]
        
        if loop_crossing != directions[1][loop_label]:
            raise ValueError("Invalid input: Must input the label's strand pd label.")

        other_potential_loop = sum(loop_crossing)/2 - (loop_label)
        if (other_potential_loop != loop_label ) and (other_potential_loop in loop_crossing) and ncrossings > 1: 
            raise ValueError("Invalid input: Cannot untwist a twisted unknot if it's not the whole link. Refer to death method instead.")

        # Determining incoming and outgoing strands + orientation + strand_type

        incoming_strand, outgoing_strand, orientation , strand_type = 0,0,0,0
        if loop_crossing[0] == loop_crossing[1] == loop_label:            # [loop_label, loop_label, new_strand_label, strand]
            orientation,strand_type = 1,1
            incoming_strand,outgoing_strand = loop_crossing[3],loop_crossing[2]
            
        elif loop_crossing[2]== loop_crossing[3] == loop_label:           # [strand, new_strand_label, loop_label, loop_label]
            orientation,strand_type = 1,-1
            incoming_strand,outgoing_strand = loop_crossing[0],loop_crossing[1]
            
        elif loop_crossing[0] == loop_crossing[3] == loop_label:          # [loop_label, strand, new_strand_label, loop_label]
            orientation,strand_type = -1,1
            incoming_strand,outgoing_strand = loop_crossing[1],loop_crossing[2]
            
        elif loop_crossing[1] == loop_crossing[2] == loop_label:          # [strand, loop_label, loop_label, new_strand_label]
            orientation,strand_type = -1,-1   
            incoming_strand,outgoing_strand = loop_crossing[0],loop_crossing[3]
            
        else:
            raise ValueError("Invalid input: Crossing is not a loop.")

        # Adjusting strands
        index_of_crossing = bad_crossings.index(loop_crossing)
        self.permute(canonical_perm(ncrossings, {index_of_crossing:ncrossings-1}), print_pd = False)
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)
        new_crossings.pop()
        if ncrossings > 1:
            other_crossing = directions[1][outgoing_strand]
            j = new_crossings.index(other_crossing)
            new_crossings[j][strand_index(new_crossings[j],outgoing_strand,1)] = incoming_strand

        self.current_min_index = 0 if not new_crossings else min([x for crossing in new_crossings for x in crossing]) 

        new_link = Link(new_crossings)
        self.links.append(new_link)

        # Computing chain map

        last_degree = self.last_degree 
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree, True, ring=self.ring) 
        self.complexes.append(final_complex)
        self.bases.append(final_bases)
            
        chain_maps = {}

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else ncrossings - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree) not in final_bases) else len(final_bases[(homological_degree, last_degree)]) 
            m = matrix(self.ring, domain_size, image_size)
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree)][image_state_index]

                    if image_state[0]== state[0][:-1]:
                    
                        remaining_circles1 = delete_circles(state, [loop_label, incoming_strand])
                        remaining_circles2 = delete_circles(image_state, [incoming_strand]) 
                        
                        if orientation == 1:
                            if state[0][-1] == 0 and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and (strand_sign(loop_label, state) == -1):
                                m[state_index,image_state_index] = 1
                        else:
                            if state[0][-1] == 1 and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])):
                                m[state_index,image_state_index] = strand_sign(loop_label,state) 

            chain_maps[homological_degree]= m.transpose()

        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)
        

        if print_pd:
            print(new_crossings)
        return new_link
 
    

    #R2

    def poke(self, strand1, strand2, parallel=-1, over=1, print_pd=True):
        
        """
        Creates a poke (Reidemeister 2 Move)
        Computes the new induced link as well as the appropriate chain maps.
        
        INPUTS: 
        strand1: First strand number in last link's pd notation
        strand2: Second strand number in last link's pd notation
        parallel : 1 if the (oriented) strands are parallel, -1 otherwise.
        over : 1 if we want strand 1 to poke on top of strand 2, -1 otherwise. 
        print_pd : If we want to output the pd code of the final link or not

        OUTPUTS: 
        Returns the last link.

        REMARK : The numbering of the strands goes as follows: outgoing parts of strand1, strand2 remain the same. We then number down mid1
        and mid2, and finally last1 and last2, in this order. For the unknot, the convention is given by (flipped labelling):
        st1,mid1,last1,st2,mid2,last2 = -1,-3,-2,-2,-4,-1
        If strand1 = strand2, the labelling is the exact same except we change strand2 for last1.
        Finally, the
        """
        
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)
        last_degree = self.last_degree

        # Adding strands 
        st1 = strand1
        st2 = strand2
        mid1 = self.current_min_index - 1
        mid2 = self.current_min_index - 2
        last1 = self.current_min_index - 3
        last2 = self.current_min_index - 4
        
        if not crossings:
            st1,mid1,last1,st2,mid2,last2 = -1,-2,-3,-3,-4,-1
            if parallel==1:
                raise ValueError("Cannot 'poke' the unknot with two parallel strands")
        else:
            if strand1==strand2:
                if parallel==1:
                    raise ValueError("Cannot 'poke' a strand into itself while being parallel")
                st2 = last1

        self.current_min_index = min(last2, mid2)
        
        # Adding crossings
        
        first_crossing = []
        second_crossing = []
        
        if over == 1 and parallel == 1 :
            first_crossing = [st2, mid1, mid2, st1]
            second_crossing = [mid2, mid1, last2, last1]
        elif over == 1 and parallel == -1:
            first_crossing = [mid2, st1, last2, mid1]
            second_crossing = [st2, last1, mid2, mid1]
        elif over == -1 and parallel == 1:
            first_crossing = [st1, st2, mid1, mid2]
            second_crossing = [mid1, last2, last1, mid2]        
        else :
            first_crossing = [st1, last2, mid1, mid2]
            second_crossing = [mid1, st2, last1, mid2] 

        new_crossings.extend([first_crossing, second_crossing])
        
        # Adjusting strands

        if crossings:
            
            directions = last_link._directions_of_edges()
            if (strand1 not in directions[0] or strand2 not in directions[0]):
                raise ValueError("Invalid input provided: One of the strands is not in the PD notation of the last link.")

            incoming_to2 = directions[1][strand2]
            j2 = crossings.index(incoming_to2)
            new_crossings[j2][strand_index(new_crossings[j2], strand2, 1)] = last2

            if strand1 != strand2:
                incoming_to1 = directions[1][strand1]
                j1 = crossings.index(incoming_to1)
                new_crossings[j1][strand_index(new_crossings[j1], strand1,1)] = last1

        new_link = Link(new_crossings)
        self.links.append(new_link)
            
            
        # Computing chain map

        last_deg = self.last_degree
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_deg, True, ring=self.ring)
        self.complexes.append(final_complex)
        self.bases.append(final_bases)

        chain_maps = {}

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else len(orientationnes) - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree) not in final_bases) else len(final_bases[(homological_degree, last_degree)]) 
            m = matrix(self.ring, domain_size, image_size)
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree)][image_state_index]

                    if image_state[0][:-2] == state[0]:
                    
                        identity_res = (0,1) if over==1 else (1,0)
                        other_res = (1,0) if over==1 else (0,1)
                        
                        remaining_circles1 = delete_circles(state, [strand1, strand2])
                        remaining_circles2 = delete_circles(image_state, [st1, st2,last1, last2, mid1, mid2])
                        
                        if image_state[0][-2:] == identity_res and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and strand_sign(strand1, state) == strand_sign(st1, image_state):
                            m[state_index, image_state_index] = 1
                        elif image_state[0][-2:] == other_res and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0]))  and strand_sign(mid1, image_state) == 1:
                            m[state_index,image_state_index] = 1

            chain_maps[homological_degree]= m.transpose()

        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)

        if print_pd:
            print(new_crossings)       
        return new_link

        

    def unpoke(self, strand1, strand2, print_pd = True):
        
        """
        Unpokes a link, opposite to the poke (R2) move. 

        INPUTS:
        strand1: The label of the first strand (the one going around the other strand). This should be the label of the strand incoming to the
                first crossing representing the poke.
        strand2: The label of the second strand (the one that remains still while the other goes around). The second strand should consist of three
                labels: two externals and one middle label. strand2 should be the label of the external strand incoming into a poke's crossing.
        print_pd : If we want to output the pd code of the final link or not.
        
        OUTPUTS:
        Returns the last link.
        """
        
        bad_crossings = self.links[-1].pd_code()
        directions = self.links[-1]._directions_of_edges()
        ncrossings = len(bad_crossings)
        
        # Recovering strands/crossings/parallel/over

        first_crossing1 = directions[1][strand1]
        mid1 = first_crossing1[((strand_index(first_crossing1, strand1, 1)+2)%4)]
        second_crossing1 = directions[1][mid1]
        last1 = second_crossing1[((strand_index(second_crossing1, mid1, 1)+2)%4)]

        first_crossing2 = directions[1][strand2]
        mid2 = first_crossing2[((strand_index(first_crossing2, strand2, 1)+2)%4)]
        second_crossing2 = directions[1][mid2]
        last2 = second_crossing2[((strand_index(second_crossing2, mid2, 1)+2)%4)]
        
        if first_crossing1 == first_crossing2 and second_crossing1 == second_crossing2 :
            parallel = 1
        elif first_crossing1 == second_crossing2 and second_crossing1 == first_crossing2:
            parallel = -1
        else:
            raise ValueError("Invalid input: The strands inputted do not represent a 'poke' move")

        over = -1 if first_crossing1[0] == strand1 else 1

        # Permuting order

        index1,index2 = bad_crossings.index(first_crossing1), bad_crossings.index(second_crossing1)

        permutation = canonical_perm(ncrossings, {index1:ncrossings-2, index2:ncrossings-1})
        self.permute(permutation, print_pd = False)

        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)[:-2]
        

        # Adjusting strands

        st1, st2 = strand1, strand2 

        if (not last1 == strand2) and (not last2==strand1): # no collapse
            
            j1 = new_crossings.index(directions[1][last1])
            j2 = new_crossings.index(directions[1][last2])
            new_crossings[j1][strand_index(new_crossings[j1], last1, 1)] = st1 
            new_crossings[j2][strand_index(new_crossings[j2], last2, 1)] = st2
            
        elif (not last1 == strand2) and (last2==strand1) : # one collapse

            j1 = new_crossings.index(directions[1][last1])
            i2 = new_crossings.index(directions[0][strand2])
            new_crossings[j1][strand_index(new_crossings[j1], last1, 1)] = st1 
            new_crossings[i2][strand_index(new_crossings[i2], strand2, -1)] = st1             

        
        elif (last1== strand2) and (not last2==strand1): # one collapse
            
            j2 = new_crossings.index(directions[1][last2])
            new_crossings[j2][strand_index(new_crossings[j2], last2, 1)] = st1

    
            
        self.current_min_index = 0 if not new_crossings else min([x for crossing in new_crossings for x in crossing]) 

        new_link = Link(new_crossings)
        self.links.append(new_link)
        
        # Computing chain map

        last_degree = self.last_degree 
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree, True, ring=self.ring) 
        self.complexes.append(final_complex)
        self.bases.append(final_bases)
            
        chain_maps = {}

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else ncrossings - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree) not in final_bases) else len(final_bases[(homological_degree, last_degree)]) 
            m = matrix(self.ring, domain_size, image_size)
            
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree)][image_state_index]
                        
                    if state[0][:-2] == image_state[0] and state[0][-1] + state[0][-2] == 1: # Appropriate states only
                    
                        remaining_circles1 = delete_circles(state, [strand1, mid2, last1])
                        remaining_circles2 = delete_circles(image_state, [st1, st2]) 
                        if not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])): # Checks that signs of unrelated circles are the same
                        
                            if ((mid1,mid2),(mid2,mid1)) in state[1] or ((mid2,mid1),(mid1,mid2)) in state[1]:
                                m[state_index,image_state_index] = -1
                            elif (((mid2, mid1),(mid1, mid2)) not in state[2]) and (((mid1, mid2),(mid2, mid1)) not in state[2]) and strand_sign(strand1, state) == strand_sign(st1, image_state):
                                m[state_index,image_state_index] = 1

            chain_maps[homological_degree]= m.transpose()

        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)

        if print_pd:
            print(new_crossings)
        return new_link

    # R3

    def slide(self, left_moving_strand, orientation_moving_strand, orientation_left_strand, orientation_right_strand, print_pd = True):
        
        """
        Does a slide move. Perspective should be a slide from TOP to BOTTOM

        INPUTS:
        left_moving_strand: Strand index of the left part of the strand that will slide.
        orientation_moving_strand: 1 if moving strand goes from left to right, -1 otherwise.
        orientation_left_strand: 1 if (bottom) left strand goes from down to up, -1 otherwise.
        orientation_right_strand: 1 if (bottom) right strand goes from down to up, -1 otherwise.
        print_pd: If we want to output the pd code of the final link or not.
        
        OUTPUTS:
        Returns the last link.
        """
        
        directions = self.links[-1]._directions_of_edges()
        bad_crossings = self.links[-1].pd_code()
        ncrossings = len(bad_crossings)
        
        # Recovering strand numbers and crossings
        # Strand notation: lms (left moving strand), mms (middle moving strand), rms (right moving strand), dls (down left strand), mls (middle left strand), uls (up left strand), drs, mrs, urs

        lms = left_moving_strand
        first_crossing = directions[ (1 if orientation_moving_strand==1 else 0) ][lms]
        mms = first_crossing[(strand_index(first_crossing, lms, orientation_moving_strand) + 2) % 4]
        second_crossing = directions[ (1 if orientation_moving_strand==1 else 0)][mms]
        rms = second_crossing[(strand_index(second_crossing, mms, orientation_moving_strand) + 2) % 4]

        over = -1 if ( (first_crossing[0] == lms and orientation_moving_strand==1) or (first_crossing[0]==mms and orientation_moving_strand==-1) ) else 1
        
        ind1 = orientation_moving_strand if over==-1 else orientation_right_strand
        ind2 = orientation_moving_strand if over==-1 else orientation_left_strand
    
        uls = first_crossing[(3 if ind1 == 1 else 1)-(1 if over==1 else 0)]  
        mls = first_crossing[(1 if ind1 == 1 else 3)-(1 if over==1 else 0)]  
        urs = second_crossing[(3 if ind2 == 1 else 1)-(1 if over==1 else 0)]  
        mrs = second_crossing[(1 if ind2 == 1 else 3)-(1 if over==1 else 0)]  

        third_crossing = directions[(1 if orientation_right_strand==-1 else 0)][mls]

        drs = third_crossing[ (strand_index(third_crossing, mls, -orientation_right_strand) + 2) % 4]
        dls = third_crossing[ (strand_index(third_crossing, mrs, -orientation_left_strand) + 2) % 4]
        
        type3 = 1 if third_crossing[1+orientation_right_strand] == mls else -1

        # Adjusting crossings

        index1, index2, index3 = bad_crossings.index(first_crossing), bad_crossings.index(second_crossing), bad_crossings.index(third_crossing) 
        new_first_crossing, new_second_crossing, new_third_crossing = [],[],[] 

        
        if over == 1:
            new_first_crossing = [dls, mms, mls, lms] if orientation_left_strand == 1 else [mls, lms, dls, mms]
            new_second_crossing = [drs, rms, mrs, mms] if orientation_right_strand == 1 else [mrs, mms, drs, rms]
        else:
            new_first_crossing = [lms, dls, mms, mls] if orientation_moving_strand == 1 else [mms, mls, lms, dls]
            new_second_crossing = [mms, drs, rms, mrs] if orientation_moving_strand == 1 else [rms, mrs, mms, drs]

        if type3 == 1:
            new_third_crossing = [mrs, urs, uls, mls] if orientation_right_strand == 1 else [ uls, mls, mrs, urs]
        else:
            new_third_crossing = [mls, mrs, urs, uls] if orientation_left_strand == 1 else [urs, uls, mls, mrs]


        permutation = {}
        final_three_crossings = []

        if over == 1 and type3 == 1: # TABLE A
            can_perm_arg = {index1:ncrossings-3,index2:ncrossings-2,index3:ncrossings-1}
            final_three_crossings = [new_second_crossing,new_first_crossing,new_third_crossing]
        elif over == -1 and type3 == 1: # TABLE B 
            can_perm_arg = {index2:ncrossings-3,index1:ncrossings-2,index3:ncrossings-1}
            final_three_crossings = [new_first_crossing,new_second_crossing,new_third_crossing]
        elif over==1 and type3 == -1: # TABLE C
            can_perm_arg = {index1:ncrossings-3,index2:ncrossings-2,index3:ncrossings-1}
            final_three_crossings = [new_second_crossing,new_first_crossing,new_third_crossing]
        else:                        # TABLE D
            can_perm_arg = {index2:ncrossings-3,index1:ncrossings-2,index3:ncrossings-1}
            final_three_crossings = [new_first_crossing,new_second_crossing,new_third_crossing]

        permutation = canonical_perm(ncrossings, can_perm_arg)
        self.permute(permutation, print_pd=False)
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)[:-3] + final_three_crossings    

        new_link = Link(new_crossings)
        self.links.append(new_link)
        
        # Computing chain map

        last_degree = self.last_degree 
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree, True, ring=self.ring) 
        self.complexes.append(final_complex)
        self.bases.append(final_bases)
            
        chain_maps = {}

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else ncrossings - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree) not in final_bases) else len(final_bases[(homological_degree, last_degree)]) 
            m = matrix(self.ring, domain_size, image_size)
            
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree)][image_state_index]
                        
                    if state[0][:-3]== image_state[0][:-3]: # Only consider appropriate states

                        remaining_circles1 = delete_circles(state, [uls, mls, dls, urs, mrs, drs, lms, mms, rms])
                        remaining_circles2 = delete_circles(image_state, [uls, mls, dls, urs, mrs, drs, lms, mms, rms])

                        if not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])): # Checks that signs of unrelated circles are the same

                            coefficient = 0

                            uls1 = strand_circle(uls, state) 
                            mls1 = strand_circle(mls, state) 
                            dls1 = strand_circle(dls, state) 
                            urs1 = strand_circle(urs, state) 
                            mrs1 = strand_circle(mrs, state) 
                            drs1 = strand_circle(drs, state) 
                            lms1 = strand_circle(lms, state) 
                            mms1 = strand_circle(mms, state) 
                            rms1 = strand_circle(rms, state) 

                            uls2 = strand_circle(uls, image_state) 
                            mls2 = strand_circle(mls, image_state) 
                            dls2 = strand_circle(dls, image_state) 
                            urs2 = strand_circle(urs, image_state) 
                            mrs2 = strand_circle(mrs, image_state) 
                            drs2 = strand_circle(drs, image_state) 
                            lms2 = strand_circle(lms, image_state) 
                            mms2 = strand_circle(mms, image_state) 
                            rms2 = strand_circle(rms, image_state) 

                            x,fx = state[0][-3:], image_state[0][-3:]
                                
                            if over == 1 and type3 == 1: # TABLE A
                                if x==fx==(0,0,0) and uls1[1] == uls2[1] and urs1[1] == urs2[1] and drs1[1]==drs2[1]:
                                    coefficient=1
                                if x==(1,0,0) and fx==(0,1,0) and uls1[1]==uls2[1] and dls1[1] == dls2[1] and drs1[1] == drs2[1]:
                                    coefficient = 1
                                if x==(1,0,0) and fx==(0,0,1) and uls1[1]==uls2[1] and dls1[1] == dls2[1] and drs1[1] == drs2[1]:
                                    coefficient = 1
                                if x==(0,1,0) and fx==(1,0,0) and uls1[1]==uls2[1] and urs1[1]==urs2[1] and dls1[1] == dls2[1]:
                                    coefficient = 1 
                                if x==fx==(1,1,0) and uls1[1]==uls2[1] and dls1[1]==dls2[1] and urs1[1]==urs2[1]:
                                    coefficient = -1
                                if x==(1,1,0) and fx==(0,1,1) and mms2[1]==1 and ( (dls1[1]==dls2[1]) or (dls1[0] in {drs1[0],urs1[0]}) ):
                                    coefficient = 1
                                if x==fx==(1,0,1) and uls1[1]==uls2[1] and rms1[1]==rms2[1] and drs1[1]==drs2[1]:
                                    coefficient=-1
                                if x==(1,0,1) and fx==(0,1,1) and mms2[1]==1 and ( (uls1[1]==uls2[1]) or (uls1[0] in {rms1[0], dls1[0]})):
                                    coefficient=-1
                                if x==(0,1,1) and fx==(1,0,1) and mms1[1] == -1 and ( (drs1[1]==drs2[1]) or (drs1[0] in {uls1[0],urs1[0]})):
                                    coefficient = 1 
                                if x==fx==(0,1,1) and mms1[1]==-1 and mms2[1] ==1:
                                    coefficient = 1
                                    
                            elif over == -1 and type3 == 1: # TABLE B 
                                if x==fx==(0,0,0) and uls1[1] == uls2[1] and urs1[1] == urs2[1] and dls1[1]==dls2[1]:
                                    coefficient = 1
                                if x==(1,0,0) and fx==(0,1,0) and uls1[1]==uls2[1] and dls1[1] == dls2[1] and drs1[1] == drs2[1]:
                                    coefficient = 1
                                if x==(1,0,0) and fx==(0,0,1) and uls1[1]==uls2[1] and dls1[1] == dls2[1] and drs1[1] == drs2[1]:
                                    coefficient = 1
                                if x==(0,1,0) and fx==(1,0,0) and uls1[1]==uls2[1] and urs1[1]==urs2[1] and dls1[1] == dls2[1]:
                                    coefficient = 1 
                                if x==fx==(1,1,0) and uls1[1]==uls2[1] and drs1[1]==drs2[1] and urs1[1]==urs2[1]:
                                    coefficient = -1
                                if x==(1,1,0) and fx==(0,1,1) and mms2[1]==1 and ( (drs1[1]==drs2[1]) or (drs1[0] in {dls1[0],uls1[0]}) ):
                                    coefficient = 1
                                if x==fx==(1,0,1) and uls1[1]==uls2[1] and rms1[1]==rms2[1] and drs1[1]==drs2[1]:
                                    coefficient=-1
                                if x==(1,0,1) and fx==(0,1,1) and mms2[1]==1 and ( (uls1[1]==uls2[1]) or (uls1[0] in {rms1[0], dls1[0]})):
                                    coefficient=-1
                                if x==(0,1,1) and fx==(1,0,1) and mms1[1] == -1 and ( (drs1[1]==drs2[1]) or (drs1[0] in {uls1[0],urs1[0]})):
                                    coefficient = 1 
                                if x==fx==(0,1,1) and mms1[1]==-1 and mms2[1] ==1:
                                    coefficient = 1                 
                                
                            elif over==1 and type3 == -1: # TABLE C
                                if x==fx==(1,0,0) and uls1[1] == uls2[1] and rms1[1] == rms2[1] and dls1[1] == dls2[1]:
                                    coefficient = 1 
                                elif x==(1,0,0) and fx==(0,1,0) and mms2[1] == 1 and ( (uls1[1] == uls2[1])   or  (uls1[0] in {rms1[0], dls1[0]}) ):
                                    coefficient = 1 
                                elif x==(0,1,0) and fx==(1,0,0) and mms1[1] == -1 and ( (dls1[1]==dls2[1]) or (dls1[0] in {uls1[0], urs1[0]}) ) :
                                    coefficient = -1 
                                elif x==fx==(0,1,0) and mms1[1]==-1 and mms2[1] == 1:
                                    coefficient = -1 
                                elif x==(0,1,0) and fx==(0,0,1) and mms1[1] == -1 and ( (uls1[1]==uls2[1]) or ( uls1[0] in {urs1[0],dls1[0]} )) :
                                    coefficient = -1  
                                elif x==fx==(0,0,1) and uls1[1]== uls2[1] and dls1[1]==dls2[1] and drs1[1] == drs2[1]:
                                    coefficient = 1 
                                elif x==(1,1,0) and fx==(1,0,1) and uls1[1]==uls2[1] and urs1[1]==urs2[1] and dls1[1] == dls2[1]:
                                    coefficient = 1 
                                elif x==(1,0,1) and fx==(0,1,1) and uls1[1]==uls2[1] and dls1[1]==dls2[1] and drs1[1]==drs2[1]:
                                    coefficient = 1
                                elif x==(0,1,1) and fx==(1,0,1) and dls1[1]==dls2[1] and uls1[1]==uls2[1] and urs1[1] == urs2[1]:
                                    coefficient=1 
                                elif x==fx==(1,1,1) and uls1[1]==uls2[1] and urs1[1]==urs2[1] and dls1[1]==dls2[1]:
                                    coefficient = -1
                            else:                        # TABLE D
                                if x==fx==(1,0,0) and uls1[1] == uls2[1] and rms1[1] == rms2[1] and dls1[1] == dls2[1]:
                                    coefficient = 1
                                elif x==(1,0,0) and fx==(0,1,0) and mms2[1] == 1 and ( (uls1[1] == uls2[1])   or  (uls1[0] in {rms1[0], dls1[0]}) ):
                                    coefficient = 1 
                                elif x==(0,1,0) and fx==(1,0,0) and mms1[1] == -1 and ( (dls1[1]==dls2[1]) or (dls1[0] in {uls1[0], urs1[0]}) ) :
                                    coefficient = -1 
                                elif x==fx==(0,1,0) and mms1[1]==-1 and mms2[1] == 1:
                                    coefficient = -1 
                                elif x==(0,1,0) and fx==(0,0,1) and mms1[1] == -1 and ( (urs1[1]==urs2[1]) or ( urs1[0] in {uls1[0],dls1[0]} )) :
                                    coefficient = -1  
                                elif x==fx==(0,0,1) and uls1[1]== uls2[1] and dls1[1]==dls2[1] and urs1[1] == urs2[1]:
                                    coefficient = 1 
                                elif x==(1,1,0) and fx==(1,0,1) and uls1[1]==uls2[1] and urs1[1]==urs2[1] and dls1[1] == dls2[1]:
                                    coefficient = 1 
                                elif x==(1,0,1) and fx==(0,1,1) and uls1[1]==uls2[1] and dls1[1]==dls2[1] and drs1[1]==drs2[1]:
                                    coefficient = 1
                                elif x==(0,1,1) and fx==(1,0,1) and dls1[1]==dls2[1] and uls1[1]==uls2[1] and urs1[1] == urs2[1]:
                                    coefficient=1 
                                elif x==fx==(1,1,1) and uls1[1]==uls2[1] and urs1[1]==urs2[1] and drs1[1]==drs2[1]:
                                    coefficient = -1                            

                            m[state_index, image_state_index] = coefficient
    
                                    
                            

            chain_maps[homological_degree]= m.transpose()
        
        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)

        if print_pd:
            print(new_crossings)

        return new_link


    # Saddle moves

    def saddle(self, strand1, strand2, print_pd=True):
        
        """
        Creates a saddle move. Computes the new link as well as induced chain maps.

        INPUTS:
        strand1 : First strand in last link's pd notation.
        strand2 : Second strand in last link's pd notation.
        print_pd : If we want to output the pd code of the final link or not
        
        OUTPUTS: 
        New last link after the saddle
        """
        
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        directions = last_link._directions_of_edges()
        if strand1 not in directions[0] or strand2 not in directions[0]:
            raise ValueError("Invalid input provided: One of the strands is not in the link")

        out_1 = directions[0][strand1]
        out_2 = directions[0][strand2]
        in_1 = directions[1][strand1]
        in_2 = directions[1][strand2]
        i1,i2 = crossings.index(out_1), crossings.index(out_2)
        j1,j2 = crossings.index(in_1), crossings.index(in_2)

        new_crossings = copy.deepcopy(crossings)

        new_crossings[j2][strand_index(new_crossings[j2], strand2,1)] = strand1 
        new_crossings[j1][strand_index(new_crossings[j1], strand1, 1)] = strand2

        new_link = Link(new_crossings)
        self.links.append(new_link)

        #Computing the induced map 

        last_degree = self.last_degree
        self.last_degree = self.last_degree - 1

        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree-1, True, ring=self.ring)

        chain_maps = {}
        
        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else len(orientationnes) - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree-1) not in final_bases) else len(final_bases[(homological_degree, last_degree-1)]) 
            m = matrix(self.ring, domain_size, image_size)
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree-1)][image_state_index]

                    if image_state[0] == state[0]:
                        remaining_circles1 = delete_circles(state, [strand1, strand2])
                        remaining_circles2 = delete_circles(image_state, [strand1, strand2])
                        
                        if not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])):
                            m[state_index,image_state_index] = 1 

            chain_maps[homological_degree]= m.transpose()

        
        self.complexes.append(final_complex)
        self.bases.append(final_bases)


        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)
        
        if print_pd:
            print(new_crossings)
        return new_link

    # Birth and death

    def birth(self, print_pd = True):
        
        """
        Births a twisted unknot.
        Virtually equivalent to the usual birth move because the homotopy equivalence is a deformation retraction.
        WARNING : Calling a birth on the unknot will first twist it.
        REMARK: The convention is to birth a twisted unknot who'se main resolution does not contain a circle + understrand. In the .twist method language, this
        corresponds to passing in (-1,-1) as the last two arguments. 

        INPUTS: 
        print_pd : If we want to output the pd code of the final link or not
        
        OUTPUTS: 
        Link object with an added twist.
        """

        was_unknot = False
        
        if not self.links[-1].pd_code():
            self.twist(-1,-1,-1, False)
            was_unknot = True 
        
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)

        # Adding strands 
        
        strand1, strand2 = self.current_min_index -1, self.current_min_index - 2 
        self.current_min_index = strand2 

        # Adding crossings

        new_crossings.append([strand1, strand2, strand2, strand1])
        new_link = Link(new_crossings)
        if was_unknot:
            self.links[-1] = new_link
        else:
            self.links.append(new_link)

        # Computing chain map
        
        last_degree = self.last_degree
        self.last_degree = last_degree + 1       
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree+1, True, ring=self.ring)
        
        if was_unknot:
            self.complexes[-1]= (final_complex)
            self.bases[-1] = (final_bases)
        else:
            self.complexes.append(final_complex)
            self.bases.append(final_bases)

        chain_maps = {}

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else len(orientationnes) - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree+1) not in final_bases) else len(final_bases[(homological_degree, last_degree+1)]) 
            m = matrix(self.ring, domain_size, image_size)
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree+1)][image_state_index]

                    if image_state[0] == (*state[0],1):
                        remaining_circles1 = delete_circles(state, [strand1, strand2])
                        remaining_circles2 = delete_circles(image_state, [strand1, strand2])
                        
                        if not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and (strand_sign(strand2, image_state) == 1):
                            m[state_index,image_state_index] = 1 

            chain_maps[homological_degree]= m.transpose()

        
        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        if was_unknot:
            self.maps[-1] = chain_morphism * self.maps[-1]
        else:
            self.maps.append(chain_morphism) 
        
        if print_pd:
            print(new_crossings)
        return new_link





    def death(self, loop_label, print_pd=True):
        
        """
        Kills a twisted unknot.
        Virtually equivalent to the usual death move because the homotopy equivalence is a deformation retraction.

        INPUTS:
        loop_label: A strand in the twisted unknot.
        print_pd : If we want to output the pd code of the final link or not
        
        OUTPUTS:
        Returns the link after the death move

        REMARK: You may only kill a twisted unknot in a link with atleast 2 components. 
        WARNING: This map may induce a minus sign. It will apply the untwist chain map on the loop label given, then the theoretical death map on the 
        remaining loop after the untwist. So be careful when birthing or killing circles to only kill the extra loop. For example, if, in your cobordism
        you birth an unknot as a twist (-2,-1) and then you do some stuff to the -1 strand before making it return to its original state, you'd want to 
        call the death map here on the -2 strand left untouched. Calling it on -1 will create an artificial negative sign.
        """
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)

        # Checking input validity

        directions = last_link._directions_of_edges()
        loop_crossing = directions[0][loop_label]
        if loop_crossing != directions[1][loop_label]:
            raise ValueError("Invalid input: Must input the label's strand pd label.")

        other_potential_loop = sum(loop_crossing)/2 - (loop_label)
        if not (other_potential_loop != loop_crossing and other_potential_loop in loop_crossing) or len(crossings) == 1: 
            raise ValueError("Invalid input: Input needs to be a twisted unknot. Link should have >= 2 crossings.")

        # Determining incoming and outgoing strands + orientation + strand_type

        other_strand, orientation , strand_type = 0,0,0
        if loop_crossing[0] == loop_crossing[1] == loop_label:            # [loop_label, loop_label, new_strand_label, strand]
            orientation,strand_type = 1,1
            other_strand = loop_crossing[3]
            
        elif loop_crossing[2]== loop_crossing[3] == loop_label:           # [strand, new_strand_label, loop_label, loop_label]
            orientation,strand_type = 1,-1
            other_strand = loop_crossing[0]
            
        elif loop_crossing[0] == loop_crossing[3] == loop_label:          # [loop_label, strand, new_strand_label, loop_label]
            orientation,strand_type = -1,1
            other_strand = loop_crossing[1]
            
        elif loop_crossing[1] == loop_crossing[2] == loop_label:          # [strand, loop_label, loop_label, new_strand_label]
            orientation,strand_type = -1,-1   
            other_strand = loop_crossing[0]
            
        else:
            raise ValueError("Invalid input: Crossing is not a loop.")

        # Adjusting strands

        index_of_crossing = crossings.index(loop_crossing)
        new_crossings.pop(index_of_crossing)
        
        self.current_min_index = min([x for crossing in new_crossings for x in crossing]) 

        new_link = Link(new_crossings)
        self.links.append(new_link)
        
        # Computing chain map

        last_degree = self.last_degree 
        self.last_degree = last_degree + 1
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree+1, True, ring=self.ring) 
        self.complexes.append(final_complex)
        self.bases.append(final_bases)
            
        chain_maps = {}

        orientationnes = last_link.orientation()
        n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
        n_plus =  0 if not crossings else len(orientationnes) - n_minus

        for homological_degree in range(-n_minus, n_plus+1):
            domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
            image_size = 0 if ((homological_degree, last_degree+1) not in final_bases) else len(final_bases[(homological_degree, last_degree+1)]) 
            m = matrix(self.ring, domain_size, image_size)
            for state_index in range(domain_size):
                for image_state_index in range(image_size):
                    state = departure_bases[(homological_degree, last_degree)][state_index]
                    image_state = final_bases[(homological_degree, last_degree+1)][image_state_index]
                    
                    resolution_copy = list(state[0])
                    resolution_copy.pop(index_of_crossing)
                    
                    if resolution_copy == list(image_state[0]):
                        remaining_circles1 = delete_circles(state, [loop_label, other_strand])
                        remaining_circles2 = delete_circles(image_state, [loop_label, other_strand]) 
                        
                        if orientation == 1:
                            if (state[0][index_of_crossing] == 0) and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and (strand_sign(loop_label, state) == -1):
                                m[state_index,image_state_index] = 1
                        else:
                            if (state[0][index_of_crossing] == 1) and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])):
                                m[state_index,image_state_index] = strand_sign(loop_label,state)


            chain_maps[homological_degree]= m.transpose()

        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)

        if print_pd:
            print(new_crossings)
        return new_link
                
        
    def permute(self, permutation, print_pd=True):
        
        """
        Computes the chain isomorphism on the canonical Khovanov chain complex induced by a permutation on the order of the crossings of the last link.
        
        INPUTS:
        permutation: A permutation sigma given by a dictionary of the form {x:y} where y = sigma(x). Only need to enter entries that change.
    
        OUTPUT:
        The new link.
        
        REMARKS: Assume that the dictionary 'permutation' is valid.
        """
        
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        ncross = len(crossings)
        invperm = {value:key for key,value in permutation.items()}
        new_crossings = [ [s for s in (crossings[index] if (index not in invperm) else crossings[invperm[index]]) ]  for index in range(ncross)]

        new_link = Link(new_crossings)

        self.links[-1] = new_link
        
        last_degree = self.last_degree
        # Computing chain map
        
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_degree, True, ring=self.ring) 

        if self.maps:
 
            departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
  
            
            chain_maps = {}
    
            orientationnes = last_link.orientation()
            n_minus = 0 if not crossings else len([x for x in orientationnes if x == -1])
            n_plus =  0 if not crossings else len(orientationnes) - n_minus
            cache = {}
            for homological_degree in range(-n_minus, n_plus+1):
                
                domain_size = 0 if ((homological_degree, last_degree) not in departure_bases) else len(departure_bases[(homological_degree, last_degree)]) 
                image_size = 0 if ((homological_degree, last_degree) not in final_bases) else len(final_bases[(homological_degree, last_degree)]) 
                
                m = matrix(self.ring, domain_size, image_size)
                
                for state_index in range(domain_size):
                    for image_state_index in range(image_size):
                        
                        state = departure_bases[(homological_degree, last_degree)][state_index]
                        image_state = final_bases[(homological_degree, last_degree)][image_state_index]

                        if state[0] == tuple(( (  ( image_state[0][i] ) if ( i not in permutation ) else (image_state[0][permutation[i]])   ) for i in range(ncross))):
                            if not (state[1].intersection(image_state[2]) or state[2].intersection(image_state[1])):
                                if state[0] in cache:
                                    m[state_index,image_state_index] = cache[state[0]]
                                else:
                                    total_factor = 1
                                    one_resolutions = []
                                    for i in range(ncross):
                                        if state[0][i]== 1:
                                            one_resolutions.append(i)
                                    num_ones = len(one_resolutions)
                                    for a_index in range(num_ones):
                                        for b_index in range(a_index+1,num_ones):
                                            a_val,b_val = one_resolutions[a_index], one_resolutions[b_index]
                                            if (a_val if a_val not in permutation else permutation[a_val]) > (b_val if b_val not in permutation else permutation[b_val]):
                                                total_factor *= -1 
                                    cache[state[0]] = total_factor 
                                    m[state_index, image_state_index] = total_factor 
                                                
                                        
                                    
    
                chain_maps[homological_degree]= m.transpose()

    
            chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        
            self.maps[-1] = chain_morphism * self.maps[-1] 

        self.complexes[-1] = final_complex 
        self.bases[-1] = final_bases
            
        if print_pd:
            print(new_crossings)
        return new_link

        
    # Chain map push

    def push(self, homological_degree, vector):

        """
        Returns the image of a vector by pushing it through the chain maps as well as the basis needed.

        INPUTS:
        homological_degree: The homological degree of the starting vector.
        vector: A vector object in SageMath.

        OUTPUTS:
        A vector representing the output through the chain maps, as well as printing the basis needed to read it.
        """

        current_vec = vector 
        
        for func in self.maps:
            current_vec = func._matrix_dictionary[homological_degree] * current_vec 

        print(f"Outputted vector: {current_vec}\n")
        print("Base in order to read output:\n---------------")
        count = 0
        for state in self.bases[-1][(homological_degree,self.last_degree)]:
            coefficient = current_vec[count] 
            if coefficient != 0:
                print(f"State:{state}\nCoefficient:{coefficient}\n---------------")

            count+=1

        return current_vec
        
            



    # Plotting 



    def next_link(self):

        """
        Returns the next link in chain based on plot_index. If at the end, circles back to the start.
        Remark: For plotting purposes, we take the absolute value of all PD strands. Hence the movie should only use negative strands.

        INPUTS:
        None

        OUTPUTS:
        Next link in the sequence in absolute value of strands.
        """
        
        actual_link = self.links[self.plot_index]
        abs_link = Link([[abs(x) for x in crossing] for crossing in actual_link.pd_code()]) # no cache, fast to compute
        self.plot_index = (self.plot_index + 1) % len(self.links)
        return abs_link


    def reset_plot_index(self):
        
        """
        Resets plot_index.

        INPUTS:
        None

        OUTPUTS:
        None
        """
        
        self.plot_index = 0      
        

# KHOVANOV HOMOLOGY HELPER FUNCTIONS

def simplified_states(link):
    
    """
    Return simplified enhanced states
    
    INPUTS:
    link : A link

    OUTPUTS:
    Enhanced Kauffman states of the link without the chords (third index in branches of circles).

    Each simplified enhanced Kauffman state is a tuple containing values corresponding to the enhanced Kauffman state as in Oleg's paper.
    If s is a state, then:
    s[0] : Is a tuple containing the crossing resolution types (tuple of n 0's and 1's where n = number of crossings in link)
    s[1] : A tuple of tuple of tuples representing negative circles.
    Each inner tuple represents a circle, each inner tuple of that represents a branch of the circle, 
    and each branch of circles looks like (strand1, strand2)
    s[2] : Same as above but for positive circles
    s[3] : Homological degree of state
    s[4] : Quantum degree of state
    """
    
    if not link.pd_code():
        return (
            ((),((),),(),0,-1),
            ((),(),((),),0,1)
        )
    states = link._enhanced_states()
    simplified= []
    for state in states:
        new_state = (state[0], tuple([tuple([tuple([x[0], x[1]]) for x in y]) for y in state[1]]),tuple([tuple([tuple([x[0], x[1]]) for x in y]) for y in state[2]]), state[3], state[4])
        simplified.append(new_state)

    return tuple(simplified)



def height_khovanov_chain_complex(link, height, base_output = False, ring=QQ):

    """
    Returns the chain complex of a certain quantum degree associated to a link following Bar Natan's convention in his paper.
    WARNING : This chain complex is only an invariant of links up to homotopy equivalence. The homology of this complex is the Khovanov homology.
    You can recover the homological dimensions from the link's Poincarré polynomial (link.khovanov_polynomial()) by the coefficient of t^x * q^{height} 
    being the dimension of H^x in the chain complex outputted by this function.

    INPUTS:
    link: A Link object.
    height: The quantum degree.
    base_output: If True, will output the base of the chain complexe's free modules. This is the same as the simplified enhanced Kauffman states ordered by 
                 degree, i.e it will output a dictionnary who'se (i,j)'th value for a key (i,j) is the states of hom degree i and qdegree j
    ring: Underlying ring.

    OUTPUTS: 
    If the base is not outputted (base_output == False), then it will return a ChainComplex object. Otherwise, it will output a tuple (chain_complex, base).
    """         
        
    crossings = link.pd_code()
    ncross = len(crossings)
    states = [(_0, set(_1), set(_2), _3, _4)
              for (_0, _1, _2, _3, _4) in simplified_states(link)]
    bases = {}  # arrange them by (i,j)
    for st in states:
        i, j = st[3], st[4]
        if j == height:
            if (i, j) in bases:
                bases[i, j].append(st)
            else:
                bases[i, j] = [st]
    complexes = {}
    for (i, j), bij in bases.items():
        if (i + 1, j) in bases:
            m = matrix(ring, len(bij), len(bases[(i + 1, j)]))
            for ii in range(m.nrows()):
                V1 = bij[ii]
                for jj in range(m.ncols()):
                    V2 = bases[(i + 1, j)][jj]
                    V20 = V2[0]
                    difs = [index for index, value in enumerate(V1[0])
                            if value != V20[index]]
                    if len(difs) == 1 and not (V2[2].intersection(V1[1]) or V2[1].intersection(V1[2])):
                         # m[ii, jj] = (-1)**sum(V1[0][x] for x in range(difs[0] + 1, ncross))
                        m[ii,jj] = (-1)**sum(V1[0][x] for x in range(difs[0]))
        else:
            m = matrix(ring, len(bij), 0)
        complexes[i] = m.transpose()
        if (i - 1, j) not in bases:
            complexes[i - 1] = matrix(ring, len(bases[(i, j)]), 0)
    if base_output:
        return (ChainComplex(complexes, base_ring=ring, check=True), bases)
    return ChainComplex(complexes, base_ring=ring,check=True)

