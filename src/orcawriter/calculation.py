import configparser
import pathlib
from jinja2 import Template

INP_TEMPLATE = r'src/orcawriter/templates/orca_input_template.txt'
SH_TEMPLATE = r'src/orcawriter/templates/orca_sh_template.txt'


class CalculationBase:
    """Base Class for a calculation object.

    Currently under development for use with the ORCA DFT package.
    However, this can be safely subclassed and used for other models and calculations.
    """

    def __init__(self, **kwargs):
        """Create an empty object, unless kwargs are passed."""
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        """Return str with the class name, and wanted attrs."""
        repr_str = f'{self.__class__.__name__}: '
        wanted_repr_attrs = ['id', 'name']
        for attr in wanted_repr_attrs:
            repr_str += f'{attr}={getattr(self, attr, __default=None)}'
        return repr_str

    def set_defaults(self, defaults_inifile: str) -> None or TypeError:
        """Set any essential parameters not passed using an .ini config file

        If the attribute is not present in the calculation object it will be set.
        Accepts .ini files, else will raise a TypeError. See Below for .ini structure


        In the base class all sections are parsed, but do not have meaning.

        .ini file structure
        [Section Heading]
        param1 = value1
        param2 = value2

        config parser makes this a dict of dicts
        {
            Heading :
                {
                    param1 : value1,
                    param2 : value2,
                    ...
                }
            Next Heading : {}
        }
        """
        defaults_inifile = pathlib.Path(defaults_inifile)
        if not defaults_inifile.suffix == '.ini':
            return TypeError('Please use an .ini file for default values.')
        config = configparser.ConfigParser()
        config.read(defaults_inifile)
        for section in config:
            for param, value in config[section].items():
                if not hasattr(self, param):
                    setattr(self, param, value)


class Calculation(CalculationBase):
    """Calculation object to represent an ORCA 4.0 calculation.

    Methods:
        to_template: renders the calculation info into a jinja2 template
        to_inp: shortcut for rendering the .inp template (using to_template)
        to_sh: shortcut for rendering the .sh template (using to_template)
    """

    def to_template(self, template: str) -> str:
        """Render a str from a jinja2 template using attr values to fill in the template.

        This str can be written to file to generate the text file desired.
        :param template: str path to a jinja2 template file
                         to render using the Calculation attrs.
        :return: output_str str
        """
        with open(template, 'r') as f:
            template_str = f.read()
        tm = Template(template_str)
        output_str = tm.render(**self.__dict__)
        return output_str

    def to_inp(self) -> str:
        """Corrects attrs and then passes INP Template automatically. Currently for ORCA 4.0

        This is merely a shortcut method, and a container for .inp specific corrections.
        Subclasses should house these corrections.
        """
        if getattr(self, 'relativistic') == 'ZORA':
            self.basis_set = f'ZORA-{self.basis_set}'
        return self.to_template(template=INP_TEMPLATE)

    def to_sh(self):
        """Corrects attrs and then passes Shell file template automatically.

        Currently for ORCA 4.0 and Yale U.
        This is merely a shortcut method, and a container for .sh specific corrections.
        Subclasses should house these corrections.
        """
        return self.to_template(template=SH_TEMPLATE)

