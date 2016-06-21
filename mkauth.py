def mkauthortex(md='F'):
    afil = [] #make affiliation list
    #affiliation 0 is Portsmouth
    afil.append('Institute of Cosmology \& Gravitation, Dennis Sciama Building, University of Portsmouth, Portsmouth, PO1 3FX, UK')
    #affiliation 1 is MPE
    afil.append('Max-Planck-Institut f\\"ur extraterrestrische Physik, Postfach 1312, Giessenbachstr., 85748 Garching, Germany')
    #affiliation 2 is LBNL
    afil.append('Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA 94720, USA')
    #affiliation 3 is CMU
    afil.append('Department of Physics, Carnegie Mellon University, 5000 Forbes Avenue, Pittsburgh, PA 15213, USA')
    #affiliation 4 is Swinburne
    afil.append('Centre for Astrophysics and Supercomputing, Swinburne University of Technology, P.O. Box 218, Hawthorn, Victoria 3122, Australia')
    #affiliation 5 is Berkeley Physics
    afil.append('Department of Physics, University of California, 366 LeConte Hall, Berkeley, CA 94720, USA')
    #affiliation 6 is Berkeley Astro
    afil.append('Department of Astronomy, University of California at Berkeley, Berkeley, CA 94720, USA')
    #affiliation 7 is Harvard
    afil.append('Harvard-Smithsonian Center for Astrophysics, 60 Garden St., Cambridge, MA 02138, USA')
    #affiliation 8 is Arizona
    afil.append('Steward Observatory, University of Arizona, 933 North Cherry Ave., Tucson, AZ 85721, USA')
    #affiliation 9 is Princeton
    afil.append('Department of Astrophysical Sciences, Princeton University, Ivy Lane, Princeton, NJ 08544, USA')
    #affiliation 10 is YCAA
    afil.append('Yale Center for Astronomy and Astrophysics, Yale University, New Haven, CT 06511, USA')
    #affiliation 11 Utah
    afil.append('Department Physics and Astronomy, University of Utah, 115 S 1400 E, Salt Lake City, UT 84112, USA')
    #affiliation 12 is Granada
    afil.append('Instituto de Astrof\\\'isica de Andaluc\\\'ia (CSIC), E-18080 Granada, Spain')
    #affiliation 13 is Case
    afil.append('Department of Astronomy, Case Western Reserve University, Cleveland, OH 44106, USA')
    #affiliation 14 is APC
    afil.append('APC, Astroparticule et Cosmologie, Universit\\\'e Paris Diderot, CNRS/IN2P3, CEA/Irfu, Observatoire de Paris, Sorbonne Paris Cit\\\'e, 10, rue Alice Domon \\& L\\\'eonie Duquet, 75205 Paris Cedex 13, France')
    #affiliation 15 is Yale physics
    afil.append('Department of Physics, Yale University, 260 Whitney Ave, New Haven, CT 06520, USA')
    #affiliation 16 CIE UAM CSIC
    afil.append('Campus of International Excellence UAM+CSIC, Cantoblanco, E-28049 Madrid, Spain')
    #affiliation 17 is IFT Madrid
    afil.append('Instituto de F\\\'isica Te\\\'orica (UAM/CSIC), Universidad Autonoma de Madrid, Cantoblanco, E-28049 Madrid, Spain')
    #affiliation 18 is CCPP NYU
    afil.append('Center for Cosmology and Particle Physics, New York University, New York, NY 10003, USA')
    #affiliation 19 is APO
    afil.append('Apache Point Observatory and New Mexico State University, P.O. Box 59, Sunspot, NM 88349, USA')
    #affiliation 20 is PSU astronomy
    afil.append('Department of Astronomy and Astrophysics, The Pennsylvania State University, University Park, PA 16802, USA')
    #affiliation 21 is is PSU IGC
    afil.append('Institute for Gravitation and the Cosmos, The Pennsylvania State University, University Park, PA 16802, USA')
    #affiliation 22 is Valencia
    afil.append('IFIC, Universidad de Valencia-CSIC, 46071, Spain')
    #affiliation 23 is Barcelona IEEC-UB
    afil.append('ICC, University of Barcelona (IEEC-UB), Marti i Franques 1, Barcelona 08028, Spain')
    #affiliation 24 is Canaries
    afil.append('Instituto de Astrof{\\\'\i}sica de Canarias (IAC), C/V{\\\'\i}a L\\\'actea, s/n, E-38200, La Laguna, Tenerife, Spain')
    #affiliation 25 is La Laguna
    afil.append('Dpto. Astrof{\\\'\i}sica, Universidad de La Laguna (ULL), E-38206 La Laguna, Tenerife, Spain')
    #affiliation 26 is BCCP
    afil.append('Berkeley Center for Cosmological Physics, LBL and Department of Physics, University of California, Berkeley, CA 94720, USA')
    #affiliation 27 is Georgia
    afil.append('National Abastumani Astrophysical Observatory, Ilia State University, 2A Kazbegi Ave., GE-1060 Tbilisi, Georgia')
    #affiliation 28 is Hubble 
    afil.append('Hubble Fellow')
    #affiliation 29 is Potsdam
    afil.append('Leibniz-Institut f\\"{u}r Astrophysik Potsdam (AIP), An der Sternwarte 16, D-14482 Potsdam, Germany')
    #affiliation 30 is LAC Paris
    afil.append('Laboratoire Astroparticule et Cosmologie, Universite Paris 7 Denis Diderot, Paris, France')
    #affiliation 31 is other Barcelona
    afil.append('ICREA \\& ICC-UB University of Barcelona, Marti i Franques 1, 08028 Barcelona, Spain')
    #affiliation 32 is Ohio State physics 
    afil.append('Department of Physics, Ohio State University, Columbus, Ohio 43210, USA')
    #affilation 33 is brazil
    afil.append('Observat\\\'orio Nacional, Rua Gal. Jos\\\'e Cristino 77, Rio de Janeiro, RJ - 20921-400, Brazil')
    #affiliation 34 is brazil lab
    afil.append('Laborat\\\'orio Interinstitucional de e-Astronomia - LineA, Rua Gal. Jos\\\'e Cristino 77, Rio de Janeiro, RJ - 20921-400, Brazil')
    #affiliation 35 is BNL
    afil.append('Brookhaven National Laboratory, Bldg 510, Upton, New York 11973, USA')
    #affiliation 36 is Paris 6
    afil.append('Institut d\'Astrophysique de Paris, Universit\\\'e Paris 6 et CNRS, 98bis Boulevard Arago, 75014 Paris, France')
    #affiliation 37 is saclay
    afil.append('CEA, Centre de Saclay, IRFU/SPP, F-91191 Gif-sur-Yvette, France')
    #affiliation 38 is APC 2
    afil.append('APC, Universit\\\'e Paris Diderot, CNRS/IN2P3, CEA/Irfu, Observatoire de Paris, Sorbonne Paris Cit\\\'e, France')
    #affiliation 39 is OSU astronomy 
    afil.append('Department of Astronomy, Ohio State University, Columbus, Ohio, USA')
    #affiliation 40 is Irvine
    afil.append('Department of Physics and Astronomy, UC Irvine, 4129 Frederick Reines Hall, Irvine, CA 92697, USA')
    #afilliation 41 is Chinese National Astronomy Observatory
    afil.append('National Astronomy Observatories, Chinese Academy of Science, Beijing, 100012, P.R. China')
    #affiliation 42 is WA
    afil.append('Department of Astronomy, University of Washington, Box 351580, Seattle, WA 98195, USA')
    #affiliation 43 is Marseille
    afil.append('Laboratoire d\'Astrophysique de Marseille, CNRS-Universite Aix-Marseille, 38 rue F. Joliot-Curie 13388 Marseille Cedex 13, France')
    #affiliation 44 is CCAPP 
    afil.append('Center for Cosmology and Astro-Particle Physics, Ohio State University, Columbus, Ohio, USA')
    #affiliation 45 is San Diego
    afil.append('Center for Astrophysics and Space Sciences, Department of Physics, University of California, 9500 Gilman Dr., San Diego, CA 92093 USA')
    #affiliation 46 is KIAS
    afil.append('Korea Institute for Advanced Study, Dongdaemun-gu, Seoul 130-722, Korea')
    #affiliation 47 is Kavli IPMU
    afil.append('Kavli Institute for the Physics and Mathematics of the Universe (WPI), The University of Tokyo Institutes for Advanced Study, The University of Tokyo, Kashiwa, Chiba 277-8583, Japan')
    #affiliation 48 is barca
    afil.append('Institut de Ci\`encies del Cosmos (ICCUB), Universitat de Barcelona (IEEC-UB), Mart\\\'i i Franqu\`es 1, E08028 Barcelona, Spain')
    #affiliation 49 is Drexel
    afil.append('Department of Physics, Drexel University, 3141 Chestnut Street, Philadelphia, PA 19104, USA ')
    #affiliation 50 is Wisconsin
    afil.append('Department of Astronomy, University of Wisconsin-Madison, 475 N. Charter Street, Madison, WI, 53706, USA')
    #affiliation 51 is Open University
    afil.append('Department of Physical Sciences, The Open University, Milton Keynes, MK7 6AA, UK')
    #affiliation 52 is Columbia
    afil.append('Department of Astronomy, Columbia University, New York, NY, 10027, USA')
    #affiliation 53 is UCL
    afil.append('University College London, Gower Street, London WC1E 6BT, UK')
    #affillation 54 is Marseille
    afil.append('CPPM, Aix-Marseille Universit\\\'e, CNRS/IN2P3, Marseille, France')
    #affiliation 55 is DFT Madrid
    afil.append('Departamento de F\\\'isica Te\\\'orica M8, Universidad Aut\\\'onoma de Madrid, E-28049 Cantoblanco, Madrid, Spain')
    #affiliation 56 is St Andrews
    afil.append('School of Physics and Astronomy, University of St Andrews, St Andrews, KY16 9SS, UK')
    #affiliation 57 is USM
    afil.append('Universit\\"ats-Sternwarte M\\"unchen, Scheinerstrasse 1, 81679 Munich, Germany')
    #affiliation 58 is  ILP (Paris)
    afil.append('Sorbonne Universit\\\'es, Institut Lagrange de Paris (ILP), 98 bis Boulevard Arago, 75014 Paris, France')
    #affiliation 59 is LPNHE Paris
    afil.append('Laboratoire de Physique Nucl\\\'eaire et de Hautes Energies, Universit\\\'e Pierre et Marie Curie, 4 Place Jussieu, 75005 Paris, France')
    #affiliation 60 is The McWilliams Center for Cosmology
    afil.append('The McWilliams Center for Cosmology, Carnegie Mellon University, 5000 Forbes Ave., Pittsburgh, PA 15213, USA')
    #affiliantion 61 is Ohio University
    afil.append('Department of Physics and Astronomy, Ohio University, 251B Clippinger Labs, Athens, OH 45701, USA')
    #affiliation 62 is La Plata
    afil.append('Facultad de Ciencias Astron\'omicas y Geof{\'\i}sicas - Universidad Nacional de La Plata. Paseo del Bosque S/N (1900). La Plata, Argentina')
    #affiliation 63 is CONICET
    afil.append('CONICET, Rivadavia 1917, 1033, Buenos Aires, Argentina')
    #affiliation 64 is NOAO
    afil.append('National Optical Astronomy Observatory, 950 N Cherry Ave, Tucson, AZ 85719, USA')
    #affiliation 65 is PITT PACC
    afil.append('PITT PACC, Department of Physics and Astronomy, University of Pittsburgh, 3941 O\'Hara Street, Pittsburgh, PA 15260, USA')
    #affiliation 66 is MPA
    afil.append('Max Planck Institut f\\"ur Astrophysik, Karl-Schwarzschild-Stra{\ss}e 1, D-85740 Garching bei M\\"unchen, Germany')
    #affiliation 67 is King's College
    afil.append('Department of Chemistry and Physics, King\'s College, 133 North River St, Wilkes Barre, PA 18711, USA')
    
    #make list for authors, form last name, first, initials, affiliation 1, affiliation 2...
    auth = []
    auth.append(('Ross', 'Ashley',' J. ', 44,0))
    auth.append(('Percival', 'Will',' J. ', 0))
    #auth.append(('Manera', 'Marc',' ', 0,53))
    auth.append(('Tojeiro', 'Rita',' ', 56))
    #auth.append(('Samushia', 'Lado',' ', 0,27))
    auth.append(('Zhao', 'Gong-Bo',' ', 41, 0))
    auth.append(('Wang', 'Yuting', ' ', 41, 0))
    auth.append(('Nichol', 'Robert',' C. ', 0))
    #auth.append(('Thomas', 'Daniel',' ', 0))
    #auth.append(('Burden', 'Angela',' ', 0))
    auth.append(('Maraston', 'Claudia',' ', 0))
    auth.append((r'S\'anchez', 'Ariel',' G. ', 1))
    auth.append(('Grieb', 'Jan',' Niklas ', 57, 1))
    auth.append(('Salazar-Albornoz', 'Salvador', ' ',57, 1))
    #auth.append(('Montesano','Francesco',' ',1))
    auth.append(('Ho', 'Shirley',' ', 3,60))
    auth.append(('Satpathy', 'Siddharth', ' ',3,60))
    auth.append(('Alam', 'Shadab', ' ',3,60))
    #auth.append(('Reid','Beth',' ',2,28,5))
    auth.append(('White','Martin',' ',2,5))
    auth.append(('Schlegel\\thanks{BOSS PI: djschlegel@lbl.gov}','David',' J. ',2))
    #auth.append(('Bailey','Stephen',' ',2))
    auth.append(('Roe','Natalie',' A. ',2))
    #auth.append(('Ross','Nicholas',' P. ',2,49))
    #auth.append(('Kazin','Eyal',' ',4))
    #auth.append(('McBride','Cameron',' K. ',7))
    auth.append(('Eisenstein','Daniel',' J. ',7))
    auth.append(('Swanson','Molly',' E. C. ',7))
    #auth.append(('Xu','Xiaoying',' ',3))
    #auth.append(('Mehta','Kushal',' T. ',8))
    #auth.append(('Skibba','Ramin',' A. ',45))
    auth.append(('Strauss','Michael',' A. ',9))
    #auth.append(('Loomis','Craig',' ',9))
    #auth.append(('Mandelbaum','Rachel',' ',3))
    #auth.append(('Gunn','James',' E. ',9))
    #auth.append(('Lupton','Robert H.',' ',3))
    auth.append(('Wake','David',' A. ',50,51))
    auth.append(('Prada','Francisco',' ',16,17,12))
    auth.append(("Rodr\\'iguez-Torres", 'Sergio', ' A. ',16,17,55)) 
    auth.append(('Zehavi','Idit',' ',13))
    #auth.append(('Guo','Hong',' ',11))
    #auth.append(('Harding','Paul',' ',13))
    auth.append(('Vargas Maga\\~na', 'Mariana',' ',14))
    #auth.append(('Hamilton','Jean-Christophe',' ',14))
    #auth.append(('Aubourg','\\\'Eric',' ',14))
    #auth.append(('Padmanabhan','Nikhil',' ',15))
    auth.append(('Cuesta','Antonio',' J. ',48))
    #auth.append(('Parejko','John',' ',15))
    #auth.append(('De Putter','Roland',' ',22,23))
    auth.append(('Scoccola','Claudia',' G. ',61,62,17,24))
    auth.append(('Seo','Hee-Jong',' ',61))
    #auth.append(('Wagner','Christian',' ',23))
    #auth.append(('Blanton','Michael',' ',18))
    auth.append(('Tinker','Jeremy',' L. ',18))
    #auth.append(('Weaver','Benjamin',' A. ',18))
    #auth.append(('Muna','Demetri',' ',18))
    #auth.append(('Mena','Olga',' ',22))
    #auth.append(('Nuza','Sebasti\\\'an',' E. ',29))
    #auth.append(('Labatie','Antione',' ',30))
    #auth.append(('Verde','Licia',' ',31))
    auth.append(('Olmstead','Matthew',' D. ',67))
    auth.append(('Dawson','Kyle',' S. ',11))
    auth.append(('Bolton','Adam',' S. ',11,64))
    auth.append(('Brownstein','Joel',' R. ',11))
    #auth.append(('Honscheid','Klaus',' ',32,44))
    auth.append(('Weinberg','David',' H. ',39,44))
    #auth.append(('Da Costa','Luiz',' N. A. ',33,34))
    auth.append(('Pan','Kaike',' ',19))
    #auth.append(('Bizyaev','Dmitry',' ',19))
    #auth.append(('Malanushenko','Viktor',' ',19))
    #auth.append(('Malanushenko','Elena',' ',19))
    #auth.append(('Oravetz','Daniel',' ',19))
    #auth.append(('Simmons','Audrey',' ',19))
    #auth.append(('Brinkmann','J.',' ',19))
    #auth.append(('Sheldon','Erin',' S. ',35))
    auth.append(('Petitjean','Patrick',' ',36))
    #auth.append(('Paris','Isabelle',' ',36))
    auth.append(('Palanque-Delabrouille','Nathalie',' ',37))
    #auth.append(('Yeche','Christophe',' ',37))
    auth.append(('Schneider','Donald',' P. ',20,21))
    #auth.append(('Kirkby','David',' ',40))
    #auth.append(('Anderson','Lauren',' ',42))
    #auth.append(('Kneib','Jean-Paul',' ',43)) 	
    #auth.append(('Howlett','Cullan',' ',0))
    auth.append(('Beutler','Florian',' ',2))
    #auth.append(('Sabiu', 'Cristiano', ' G. ',46))
    auth.append(('Chuang', 'Chia-Hsun', ' ',17, 29))
    auth.append(('Saito', 'Shun', ' ', 47, 66)) 
    #auth.append(('Bhardwaj','Vaishali',' ',42,2))
    auth.append(('Price-Whelan','Adrian',' M. ',9))
    #auth.append(('Escoffier','Stephanie',' ',54))
    auth.append(('Pellejero-Ibanez', 'Marcos', ' ', 24, 25))
    auth.append(('Slosar', 'An\\v{z}e',' ', 35))
    auth.append(('Vazquez', 'Jose Alberto',' ', 35))
    auth.append(("Gil-Mar\\'in", "H\\'ector" , ' ', 58 , 59))
    auth.append(('Ata', 'Metin', ' ', 29))
    auth.append(('Zhai', 'Zhongxu', ' ', 18))
    auth.append(('Kitaura', 'Francisco', ' ', 29, 2, 6))
    auth.append(('Comparat', 'Johan', ' ', 16, 17, 55))
    auth.append(('Wood-Vasey', 'W. Michael', ' ', 65))
    auth.append(('Streblyanska', 'Alina', ' ', 24,25))
    auth.append(('Bahcall', 'Neta', ' A. ', 9))
    auth.append(("L\\'opez-Corredoira", "Mart\\'in", ' ', 24,25)) 
    auth.append(("Rubi\\~{n}o-Mart\\'in", "Jose Alberto", ' ', 24,25))
    
    auth.sort()
    afn = []
    for i in range(3,len(auth[0])):
        afn.append(auth[0][i])
    for i in range(1,len(auth)):
        for j in range(3,len(auth[i])):
            aff = auth[i][j]
            s = 0
            for k in range(0,len(afn)):
                if afn[k] == aff:
                    s = 1
            if s == 0:
                afn.append(aff)
    fo = open('authlist.tex','w')
    fo.write('\\author['+auth[0][1][0]+'. '+ auth[0][0]+' et al.]')
    fo.write('{\parbox{\\textwidth}{\n')
    for i in range(0,len(auth)):
        affal = []
        for k in range(3,len(auth[i])):
        	aff = auth[i][k]
        	for j in range(0,len(afn)):
        		if afn[j] == aff:
        			ind = j+1
        	affal.append(ind)		
#         aff0 = auth[i][3]
#         for j in range(0,len(afn)):
#             if afn[j] == aff0:
#                 ind = j+1
#         affal.append(ind)		
#         if len(auth[i]) > 4:
#             aff1 = auth[i][4]
#             for j in range(0,len(afn)):
#                 if afn[j] == aff1:
#                     ind1 = j+1
#             affal.append(ind1)		
#             #fo.write(','+str(ind1))
#         if len(auth[i]) > 5:
#             aff2 = auth[i][5]
#             for j in range(0,len(afn)):
#                 if afn[j] == aff2:
#                     ind2 = j+1
#             affal.append(ind2)
            #fo.write(','+str(ind2))
        #affal.sort()
        if md == 'F':		
            fo.write(auth[i][1]+auth[i][2]+auth[i][0])#+'$^{'+str(affal[0]))
        if md == 'i':
            fo.write(auth[i][1][0]+'.'+auth[i][2]+auth[i][0])#+'$^{'+str(affal[0]))
        
        fo.write('$^{')
        for k in range(0,len(affal)):
            if k == len(affal)-1:
                fo.write(str(affal[k]))
            else:
                fo.write(str(affal[k])+',')
        if i != len(auth)-1:
            fo.write('}$,\n')
        else:
            fo.write('}$\n')	
    fo.write(' } \\vspace*{4pt} \\\ \n')
    print afn
    for i in range(0,len(afn)):
        fo.write('\\scriptsize $^{'+str(i+1)+'}$ '+afil[afn[i]]+'\\vspace*{-2pt} \\\ \n')
    fo.write('}')
    fo.close()


if __name__ == '__main__':

    mkauthortex()
