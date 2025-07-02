# Assumes SageMath is implemented. Suggested use on Jupyter.


def height_khovanov_chain_complex(link, height, ring=QQ):
  """
  Returns the associated chain complex to a link at a certain quantum degree. 

  INPUTS:
  - link : The Link object in Sagemath.
  - height : Quantum degree of the Chain Complex.
  - ring : Ring on which the computation is done (May not work if ring is not a field of characteristic 0)

  OUTPUT:
  - A ChainComplex object representing the usual Bar-Natan's construction for a complex associated to a link at a specific quantum degree.
  The ChainComplex is represented as a collection of matrices.

  """
        crossings = link.pd_code()
        ncross = len(crossings)
        states = [(_0, set(_1), set(_2), _3, _4)
                  for (_0, _1, _2, _3, _4) in link._enhanced_states()]
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
        return ChainComplex(complexes)
