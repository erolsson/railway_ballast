from __future__ import print_function
import os
try:
    import distro
except ImportError:
    import platform as distro

package_path = os.path.dirname(__file__)

if distro.linux_distribution()[0] == 'Ubuntu':
    abq = 'singularity exec --nv ' + os.path.expanduser('~/imgs/sing/abaqus-2018-centos-7.img') + \
          ' vglrun /opt/abaqus/2018/Commands/abq2018'
else:
    abq = '/scratch/users/erik/SIMULIA/CAE/2018/linux_a64/code/bin/ABQLauncher'