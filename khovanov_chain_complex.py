import copy


# HELPER FUNCTIONS

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

    raise ValueError("get your crossings right or smth")


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

    if not len(state[0]):
        return -1 if len(state[1]) else 1

    for circle in state[1]:
        for branch in circle:
            if strand in branch:
                return -1

    return 1

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
    -untwist : Undo an R1 move (TODO)
    -poke : Does an R2 move
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

    Remark: ONLY use negative strands in PD notation.
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
        
        self.links = [link] 
        self.current_min_index = 0 if not link.pd_code() else min([x for crossing in link.pd_code() for x in crossing]) #
        self.plot_index = 0
        self.last_degree = starting_qdeg 
        self.maps = [] 
        self.ring = ring 
        self.bases = [] #
        self.complexes = [] 

        initial_complex, initial_base = height_khovanov_chain_complex(link, starting_qdeg, True, ring=self.ring)
        self.bases.append(initial_base)
        self.complexes.append(initial_complex)

    # REIDEMEISTER MOVES

    # R1 
    
    def twist(self, strand, orientation = 1, strand_type = 1):
        
        """
        Creates a twist (Reidemeister 1 Move)
        Computes the new induced link as well as the appropriate chain maps.
        
        INPUTS: 
        strand: Strand number of the last link in PD notation
        orientation : If we want the twist's loop's 0-resolution to contain a loop or not. If orientation = 1, it will, otherwist it won't. 1 by default.
        strand_type : If we want the loop's crossing to be overstrand (strand_type = 1) or understrand (strand_type = -1)

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

        # GENERAL CASE

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

        # Computing induced chain map
        
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

                    remaining_circles1 = delete_circles(state, [strand])
                    remaining_circles2 = delete_circles(image_state, [strand, loop_label, new_strand_label])
                    
                    if orientation == 1:
                        if image_state[0] == (*state[0],0) and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])):
                            m[state_index,image_state_index] = -1 if (((loop_label, loop_label),) in image_state[2]) else 1
                    else:
                        if image_state[0] == (*state[0],1) and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and (((loop_label, loop_label),) in image_state[2]):
                            m[state_index,image_state_index] = 1 

            chain_maps[homological_degree]= m.transpose()

            
            
        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)

        
        return new_link

    # def untwist(self, loop_label):
    #     # last_link = self.links[-1]
    #     # crossings = last_link.pd_code()
    #     # new_crossings = copy.deepcopy(crossings)

    #     # directions = last_link._directions_of_edges()
    #     # loop_crossing = directions[0][loop_label]
    #     # if loop_crossing != directions[1][loop_label]:
    #     #     raise ValueError("Invalid input : Must input the label's strand pd index")

    #     # if (sum(loop_crossing)/2 - *(loop_label)) in loop_crossing: 

     
    #     return None

    

    #R2

    def poke(self, strand1, strand2, parallel=-1, over=1):
        
        """
        Creates a poke (Reidemeister 2 Move)
        Computes the new induced link as well as the appropriate chain maps.
        
        INPUTS: 
        strand1: First strand number in last link's pd notation
        strand2: Second strand number in last link's pd notation
        parallel : 1 if the (oriented) strands are parallel, -1 otherwise.
        over : 1 if we want strand 1 to poke on top of strand 2, -1 otherwise. Only relevant if strand1==strand2 otherwise switching order is enough

        OUTPUTS: 
        Returns the last link.

        REMARK : The numbering of the strands goes as follows: outgoing parts of strand1, strand2 remain the same. We then number down mid1
        and mid2, and finally last1 and last2, in this order. Should there be some collapses, we take the minimum label. The labelling is flipped
        for the unknot, so we suggest looking up the PD notation to start with that in case your movie starts with the unknot. 
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
            st1,mid1,last1,st2,mid2,last2 = -1,-3,-2,-2,-4,-1
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

        print(new_crossings)
        new_link = Link(new_crossings)
        self.links.append(new_link)
            
            
        # Computing chain map

        last_deg = self.last_degree
        departure_complex, departure_bases = self.complexes[-1], self.bases[-1]
        final_complex, final_bases = height_khovanov_chain_complex(new_link, last_deg, True, self.ring)
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
                    
                    identity_res = (*state[0],0,1) if over==1 else (*state[0],1,0)
                    other_res = (*state[0],1,0) if over==1 else (*state[0],0,1)
                    
                    remaining_circles1 = delete_circles(state, [strand1, strand2])
                    remaining_circles2 = delete_circles(image_state, [st1, st2,last1, last2, mid1, mid2])
                    
                    if image_state[0] == identity_res and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])) and strand_sign(strand1, state) == strand_sign(st1, image_state):
                        m[state_index, image_state_index] = 1
                    elif image_state[0] == other_res and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0]))  and strand_sign(mid1, image_state) == 1:
                        m[state_index,image_state_index] = 1

            chain_maps[homological_degree]= m.transpose()

        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)


        
        return new_link

    # def unpoke(self):
    #     """
    #     TODO
    #     """
    #     return None



    # Saddle moves

    def saddle(self, strand1, strand2):
        
        """
        Creates a saddle move. Computes the new link as well as induced chain maps.

        INPUTS:
        strand1 : First strand in last link's pd notation.
        strand2 : Second strand in last link's pd notation.

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

                    remaining_circles1 = delete_circles(state, [strand1, strand2])
                    remaining_circles2 = delete_circles(image_state, [strand1, strand2])
                    
                    if image_state[0] == state[0] and not (remaining_circles1[0].intersection(remaining_circles2[1]) or remaining_circles1[1].intersection(remaining_circles2[0])):
                        m[state_index,image_state_index] = 1 

            chain_maps[homological_degree]= m.transpose()

        
        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)
        self.complexes.append(final_complex)
        self.bases.append(final_bases)
        
        return new_link

    # Birth and death

    def twist_birth(self):
        
        """
        Births a twisted unknot.
        Virtually equivalent to the usual birth move because the homotopy equivalence is a deformation retraction.
        WARNING: Calling a twist_birth on the unknot will simply untwist it.
        """
        return None





    def twist_death(self):
        
        """
        Kills a twisted unknot.
        Virtually equivalent to the usual death move because the homotopy equivalence is a deformation retraction.
        """
        return None




    # Plotting 



    def next_link(self):

        """
        Returns the next link in chain based on plot_index. If at the end, circles back to the start.
        Remark: For plotting purposes, we take the absolute value of all PD strands. Hence the movie should only use negative strands.
        """
        
        actual_link = self.links[self.plot_index]
        abs_link = Link([[abs(x) for x in crossing] for crossing in actual_link.pd_code()]) # no cache, fast to compute
        self.plot_index = (self.plot_index + 1) % len(self.links)
        return abs_link


    def reset_plot_index(self):
        """
        Resets plot_index.
        """
        self.plot_index = 0
                
        


# KHOVANOV HOMOLOGY 


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
                        m[ii, jj] = (-1)**sum(V2[0][x] for x in range(0,difs[0]))
                        # Here we have the matrix constructed, now we have to put it in the dictionary of complexes
        else:
            m = matrix(ring, len(bij), 0)
        complexes[i] = m.transpose()
        if (i - 1, j) not in bases:
            complexes[i - 1] = matrix(ring, len(bases[(i, j)]), 0)
    if base_output:
        return (ChainComplex(complexes), bases)
    return ChainComplex(complexes)






