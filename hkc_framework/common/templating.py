from string import Template

def templating(self, tmpfilePath, writePath, fileName, parameters):
    """
    Purpose: Handles the creation and templating of files.

    Args:
        tmpfilePath = Path to tmp file that contains the necessary styling for templating.
        writePath = Write path for templated file.
        fileName = name of file and extension.
        parameters = substuting dictionary.
    """

    filein = open(tmpfilePath, "r")
    src = Template(filein.read())
    subbedtemplate = src.substitute(parameters)
    writeFile = open(writePath + "/" + fileName, "w")
    writeFile.write(subbedtemplate)
    writeFile.close()




