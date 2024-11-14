import numpy as np
from numpy.lib.recfunctions import append_fields
from list_isos import isos_from_net

def blend_comps(comp_surf, configs, multiplier=1e-6):
    """Create a structured array of blended compositions.
    
    Parameters
    ----------
    comp_surf : np.array
        The surface composition.
    configs : list of tuples
        Each tuple contains the transition fractional external mass coordinate,
        the width of the transition in mass coordinate space, and a one-element
        structured array of the composition to blend to through this transition.
        The list should be ordered by increasing transition external mass
        coordinate (ending with the core composition).
    multiplier : float, optional
        Relative steepness of transitions. At each transition, the difference
        in external fractional mass coordinate between the two composisitions is
        computed as this value multiplied by the external fractional mass at the
        transition. The default is 1e-6.
    
    Returns
    -------
    np.array
        A structured array with the same field names as the input compositions
        and with the first field being the fractional external mass coordinate 
        (xq). The fractional external mass coordinate will be 0 for the surface
        composition and 1 for the core composition, which should be the last
        composition in the list.
    """
    ###############################
    # Set up the structured array #
    ###############################
    dt = np.dtype([('xq', float)] + [(name, float) for name in comp_surf.dtype.names])
    # each composition gets a starting and a stopping point, so two points
    # per composition, which is the length of the configs list plus one more
    # for the surface composition.
    length = 2 * len(configs) + 2
    data = np.zeros(length, dtype=dt)

    def xq_comp_tuple(xq, comp):
        """Return a tuple of the external mass and the composition values.
        
        Main purpose is to create rows of the main structured array, which
        have a external fractional mass coordinate followed by the isotopic
        mass fractions.
        
        Parameters
        ----------
        xq : float
            The fractional external mass coordinate.
        comp : np.array
            The composition.
        
        Returns
        -------
        tuple
        """
        return tuple([xq] + [comp[name] for name in comp.dtype.names])
    
    # set core and surface compositions
    data[0] = xq_comp_tuple(0, comp_surf)
    data[-1] = xq_comp_tuple(1, configs[-1][-1])

    # iterate through the transition points and blend between compositions
    prev_comp = comp_surf
    for i, (xq0, dxq, comp) in enumerate(configs):
        # each transition gets two points. The first is the outer
        # composition and the second is the inner composition. Somewhat
        # arbitrarily choose the width of the transition to be 1e-6 the
        # fractional external mass coordinate at the transition to make
        # sharp, but not discontinuous transitions.

        dxq = multiplier * xq0
        # old composition
        prev_comp = comp_surf if i == 0 else configs[i-1][-1]
        data[2 * i + 1] = xq_comp_tuple(xq0-dxq, prev_comp)
        # new composition
        if i <= len(configs) - 1:
            data[2 * i + 2] = xq_comp_tuple(xq0+dxq, comp)

    return data

def alternate_iso(iso, isotopes, comp_array, use_max=True):
    """Find an alternate isotope in a composition array.
    
    If an isotope is not in the composition array, find an isotope of the same element that is in the array and has the largest total mass and return it
    as a string. If no suitable replacement is found, raise an error if `use_max` is `False` or return the most common isotope of the element with the highest atomic number in the model.
    
    Parameters
    ----------
    iso : str
        The isotope to find an alternate for.
    isotopes : list of str
        The list of acceptable isotopes.
    comp_array : np.array
        The structured array of compositions. The first field must be the
        fractional external mass coordinate (xq). The other fields must be
        isotopic mass fractions. This is the kind of array that is returned by
        `blend_comps`.
    use_max : bool, optional
        If `True`, use the element with the largest atomic number and the
        isotope of that element with the largest total mass in the model as a
        backup. If `False`, raise an error if no suitable replacement is found.

    Returns
    -------
    str
        The name of the alternate isotope.

    Raises
    ------
    ValueError
        If `use_max` is `False` and no suitable replacement is found.
    """
    element = "".join([c for c in iso if c.isalpha()])
    alternates = []
    for name in comp_array.dtype.names:
        if element == "".join([c for c in name if c.isalpha()]) and name in isotopes:
            alternates.append(name)
    if alternates:
        # find the isotope with the largest total mass in the model
        dqs = np.diff(comp_array['xq'])
        total_masses = [np.sum(comp_array[alt][:-1] * dqs) for alt in alternates]
        max_iso = alternates[np.argmax(total_masses)]
        return max_iso
    elif use_max:
        # if no suitable replacement is found, use the element with the largest 
        # atomic number and the isotope of that element with the largest total
        # mass in the model
        max_element = isotopes[-1].strip("1234567890")
        # find all isotopes of the element with the largest atomic number in
        # the model (there may be none)
        alternates = []
        for name in comp_array.dtype.names:
            if max_element == "".join([c for c in name if c.isalpha()]):
                alternates.append(name)
        if len(alternates) == 0:
            return isotopes[-1]
        elif len(alternates) == 1:
            return alternates[0]
        # find the isotope with the largest total mass in the model
        else:
            dqs = np.diff(comp_array['xq'])
            total_masses = [np.sum(comp_array[alt][:-1] * dqs) for alt in alternates]
            return alternates[np.argmax(total_masses)]
    else:
        raise ValueError(f"Isotope {iso} not in composition array and there are no suitable replacements.")

def make_composition_file(comp_array, net, filename):
    """Write a composition array to a file.

    Must have a structured array of xq coordinates and isotopic mass fractions,
    like what is created by `blend_comps`.

    Parameters
    ----------
    comp_array : np.array
        The structured array of compositions. The first field must be the
        fractional external mass coordinate (xq). The other fields must be
        isotopic mass fractions. This is the kind of array that is returned by
        `blend_comps`.
    net : str
        The name of the net file to use. This is used to get the list of isotopes in the composition array. Note that all isotopes in the net file must be present in the composition array.
    filename : str
        The name of the file to write the composition to. If the file already exists, it will be overwritten.
    
    Returns
    -------
    None
    """
    isotopes = isos_from_net(net)
    for iso in comp_array.dtype.names[1:]:
        if iso not in isotopes:
            alt = alternate_iso(iso, isotopes, comp_array, use_max=True)
            # dump this isotope's mass into the max isotope
            print(f"Isotope {iso} not in composition array. Dumping into {alt}.")
            if alt in comp_array.dtype.names:
                comp_array[alt] += comp_array[iso]
            else:
                # update comp_array to have new column with the alternate isotope
                comp_array = append_fields(comp_array, alt, comp_array[iso], usemask=False)
    header = f"{len(comp_array)} {len(isotopes)}\n"
    contents = header
    for row in comp_array:
        contents += f"{row['xq']:.15e}"
        for iso in isotopes:
            # if isotope is in composition array, add it to the line
            # if not, add 0
            if iso in comp_array.dtype.names:
                contents += f" {row[iso]:.8e}"
            else:
                contents += " 0.00000000e+00"
        contents += "\n"
    with open(filename, 'w') as f:
        f.write(contents.rstrip())
