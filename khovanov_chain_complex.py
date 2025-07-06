import copy


# HELPER FUNCTIONS

def strand_index(crossing, strand, direction):
    
    """
    Given a crossing in PD notation and a strand in that crossing, it will
    return the index (0,1,2 or 3) of the strand in the crossing. 
    This assumes the strand is in the crossing.
    If direction == 1, it will return the index at which the strand is incioming.
    Else, it will return the index at which the strand is outgoing.
    
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


# MOVIE CLASS

class Movie():
    
    """
    Class representing movie cobordisms. Explanation of attributes below:

    
    -links : The links of the movie in order as a list.
    -current_min_index : The minimum index used in the PD notations for all links
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
        self.current_min_index = 0 if not link.pd_code() else min([x for crossing in link.pd_code() for x in crossing])
        self.plot_index = 0
        self.last_degree = starting_qdeg
        self.maps = []
        self.ring = ring
        self.bases = []
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
            
            out_going_from = directions[0][strand]
            incoming_to = directions[1][strand]
            i,j = crossings.index(out_going_from), crossings.index(incoming_to)
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
                    if orientation == 1:
                        if image_state[0] == (*state[0],0) and not (state[1].intersection(image_state[2]) or state[2].intersection(image_state[1])):
                            m[state_index,image_state_index] = -1 if (((loop_label, loop_label),) in image_state[2]) else 1
                    else:
                        if image_state[0] == (*state[0],1) and not (state[1].intersection(image_state[2]) or state[2].intersection(image_state[1])) and (((loop_label, loop_label),) in image_state[2]):
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

    def poke(self, strand1, strand2):
        """
        TODO
        """
        last_link = self.links[-1]
        crossings = last_link.pd_code()
        new_crossings = copy.deepcopy(crossings)
        last_degree = self.last_degree

        # Adding crossings

    
        # Adjusting strands


        # Computing chain map

        

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
                    if image_state[0] == state[0] and not (state[1].intersection(image_state[2]) or state[2].intersection(image_state[1])):
                        m[state_index,image_state_index] = 1 

            chain_maps[homological_degree]= m.transpose()

        
        chain_morphism = Hom(departure_complex,final_complex)(chain_maps)
        self.maps.append(chain_morphism)
        self.complexes.append(final_complex)
        self.bases.append(final_bases)
        
        return new_link
        



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

        # if not link.pd_code():
        #     if height**2 == 1:
        #         return ChainComplex({0: matrix(ring, 0,1)})
        #     else:
        #         return ChainComplex({0 : matrix(ring, 0,0)})
                
                
        
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
