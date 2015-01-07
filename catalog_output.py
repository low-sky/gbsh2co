import astropy.io.fits as fits
from astropy.table import Table,Column
from astropy.coordinates import SkyCoord
import astropy.units as u



t=Table.read('s2_450.fits',format='fits')
t.add_column(Column(data=(['SCIENCE']*len(t)),name='type'))
t.add_column(Column(name='Name',data=(['                       ']*len(t))))
t.add_column(Column(name='RAString',data=(['                       ']*len(t))))
t.add_column(Column(name='DECString',data=(['                       ']*len(t))))

for idx,source in enumerate(t):
    c = SkyCoord(ra=source['RA']*u.degree,
                 dec=source['DEC']*u.degree,
                 frame='icrs')
    t['RAString'][idx] = c.ra.to_string(sep=':',pad='True',alwayssign=True,precision=2,unit=u.hourangle)
    t['DECString'][idx] = c.dec.to_string(sep=':',pad='True',alwayssign=True,precision=2)
    
    strval = c.galactic.to_string('decimal')
    if c.galactic.b.value<0:
        t['Name'][idx] = 'G{0}{1}'.format(*strval.split(' '))
    else:
        t['Name'][idx] = 'G{0}+{1}'.format(*strval.split(' '))

t['Type'] = 'SCIENCE'
t2 = t['Type','Name','RAString','DECString']
t2.write('h2cocat.txt',format='ascii')
