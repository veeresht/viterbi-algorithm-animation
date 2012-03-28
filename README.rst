The Viterbi Algorithm - Illustrated!
====================================

This software enables the generation of illustrations for the Viterbi Algorithm 
decoding of convolutional codes using Python.

A demo of the illustration created using this software can be found here_.

.. _here: http://veeresht.info/blog/viterbi-algorithm-illustrated/

Requirements
------------
- Python 2.7 or above
- NumPy 1.6 or above
- SciPy 0.10 or above
- Matplotlib 1.1 or above
- Cython 0.15 or above
- CommPy 0.1

Note: There are some modifications that need to be done in Matplotlib's **animation.py** file. 
In subclass::

    class ArtistAnimation(TimedAnimation)

comment the function::

    def _pre_draw(self, framedata, blit) 

as follows::
    
    def _pre_draw(self, framedata, blit):
        '''
        Clears artists from the last frame.
        '''
        """
        if blit:
            # Let blit handle clearing
            self._blit_clear(self._drawn_artists, self._blit_cache)
        else:
            # Otherwise, make all the artists from the previous frame invisible
            for artist in self._drawn_artists:
                artist.set_visible(False)
        """
