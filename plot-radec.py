import fitsio
import matplotlib.pyplot as mplt

data = fitsio.read('psf-data.fits.gz')

fig, ax = mplt.subplots()

ax.set(xlabel='RA', ylabel='DEC')
ax.scatter(
    data['ra'],
    data['dec'],
    s=0.1,
)

# fig.savefig('radec.png')
mplt.show()
