      module m_lmlabel
      integer,parameter:: lmlabel_lmax=2
      character(5):: lmlabel1(9)
         data lmlabel1 /'s','px','pz','pz','dxy','dyz','dz2','dxz','dx2y2'/

      character(5):: lmlabel2(-lmlabel_lmax:lmlabel_lmax,0:lmlabel_lmax) 
     
        data lmlabel2/
     .   '','','s','','',
     .   '','py','pz','px','',
     .   'dxy', 'dyz', 'dz2', 'dxz', 'dx2y2'/

      end module m_lmlabel



