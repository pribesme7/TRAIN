from parseTFS import parseTFS
import numpy as np
import matplotlib.pyplot as plt


def get_f_terms_for(R11,R12,R21,R22,betx,bety,alfx,alfy):
    NORMR = np.sqrt(1.0 + R11 * R22- R12 * R21);

    R11 /= NORMR;
    R12 /=  NORMR;
    R21 /=  NORMR;
    R22 /=  NORMR;

    #!--- Gb is actually inv(Gb)
    Ga11 = 1.0 / np.sqrt(betx);
    Ga22 = np.sqrt(betx);
    Ga21 = alfx / np.sqrt(betx);
    Ga12 = 0.0;
    Gb21 = -alfy / np.sqrt(bety);
    Gb12 = 0.0;
    Gb11 = np.sqrt(bety);
    Gb22 = 1.0 / np.sqrt(bety);

    CP11 =  R22 * Gb11 - R12 * Gb21;
    CP12 =  R22 * Gb12 - R12 * Gb22;
    CP21 =  - R21 * Gb11 + R11 * Gb21;
    CP22 =  - R21 * Gb12 + R11 * Gb22;

    C11 = Ga11 * CP11 + Ga12 * CP21;
    C12 = Ga11 * CP12 + Ga12 * CP22;
    C21 = Ga21 * CP11 + Ga22 * CP21;
    C22 = Ga21 * CP12 + Ga22 * CP22;

    GAMMA = np.sqrt(1.0 - (C11 * C22 - C12 * C21));

    elem_name_f1001r = ( C12 - C21)/4.0/GAMMA;# !--- F1001R
    elem_name_f1001i = ( C11 + C22)/4.0/GAMMA;# !--- F1001I
    elem_name_f1010r = (-C12 - C21)/4.0/GAMMA;# !--- F1010R
    elem_name_f1010i = ( C11 - C22)/4.0/GAMMA;# !--- F1010I

    return elem_name_f1001r+1j*elem_name_f1001i,elem_name_f1010r+1j*elem_name_f1010i
    


if __name__=="__main__":
    plt.figure(10)

    data = parseTFS('MAD_PART/MDcycle/train_rmatr.optb');
    for key in data.keys():
        data[key] = np.array(data[key])

    #F1001,F1010 = get_f_terms_for(data['R11'],data['R12'],data['R21'],data['R22'],data['BETX'],data['BETY'],data['ALFX'],data['ALFY'])

    #plt.plot(data['S'],data['R11'],'--')
    #plt.plot(data['S'],data['R12'],'--')
    #plt.plot(data['S'],data['R21'],'--')
    #plt.plot(data['S'],data['R22'],'--')

    #data = parseTFS('b1.twiss');
    #for key in data.keys():
    #    data[key] = np.array(data[key])

    F1001,F1010 = get_f_terms_for(data['R11'],data['R12'],data['R21'],data['R22'],data['BETX'],data['BETY'],data['ALFX'],data['ALFY'])

    plt.plot(data['S'],np.real(F1001),'-r')
    plt.plot(data['S'],np.imag(F1001),'--g')
    plt.plot(data['S'],np.real(F1010),'-',color='purple')
    plt.plot(data['S'],np.imag(F1010),'--b')


    plt.show()

