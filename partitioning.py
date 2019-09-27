from math import pi



class Partitioning:

    def __init__(self, conc_184w_silicate, D, mass_silicate, mass_metal):

        self.conc_184w_silicate = float(conc_184w_silicate)
        self.D = float(D)
        self.mass_silicate = float(mass_silicate)
        self.mass_metal = float(mass_metal)


    def __calc_mass_transfer(self):

        numerator = self.D * (self.conc_184w_silicate)
        denominator = (1 / self.mass_metal) + (self.D * (1 / self.mass_silicate))

        m = numerator / denominator

        return m

    def __verify_D(self):

        m = self.__calc_mass_transfer()

        numerator = m / self.mass_metal
        denominator = self.conc_184w_silicate - (m / self.mass_silicate)

        D = numerator / denominator

        return D

    def partition(self, metal):

        m = self.__calc_mass_transfer()

        metal.conc_184w_metal = (m / self.mass_metal)
        self.conc_184w_silicate = self.conc_184w_silicate - (m / self.mass_silicate)

        return metal.conc_184w_metal, self.conc_184w_silicate


class Metal:

    def __init__(self, density, radius, conc_184w_metal=0):

        self.conc_184w_metal = conc_184w_metal
        self.volume = (4/3) * pi * (radius)**3
        self.mass = density * self.volume
