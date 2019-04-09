from astropy import units as u
import numpy as np

def getRCont(h, rhoCrust, rhoMantle, outputUnits=None):
    rCont = (h*rhoCrust)/(rhoMantle - rhoCrust)
    if outputUnits is not None:
        rCont = rCont.to(outputUnits)

    return rCont

headerDefault ="""
\\begin{table}[H]
    \centering
    \caption{Caption}
    \label{tab:my_label}
    \\begin{tabular}{}
"""
footerDefault = """
    \end{tabular}
\end{table}
"""

def npArray2LatexTable(array, savename, header=headerDefault, footer=footerDefault, headerComment=False):
    """
    Simple function to convert numpy array to a latex table component
    :param array: Array to convert
    :param savename: Name to give to generated txt file
    :param header: Header string to give file. Default creates the beginning of a table
    :param footer: Footer string to give file. Default creates the end of a table
    :param headerComment: Boolean to say if header should be commented or not
    :return None:
    """
    if headerComment is False:
        comments=""
    else:
        comments="#"

    np.savetxt(savename, array, fmt="%s", delimiter="\t&\t", newline="      \\\ \n\hline \n",
               header=header, footer=footer, comments=comments)
