"""
In the contex of multi-fish sequential fish, different genes are targeted with the same barcode. 
Although different barcode are mixed in such a way it is possible to decode which is which using the information found on several cycles.  

    **Exemple**
    BarcodeA -> Gene1, Gene2, Gene3
    BarcodeB -> Gene3, Gene4
    BarcodeC -> Gene1, Gene4

    Gene1 = A and C
    Gene2 = A only
    Gene3 = A and B 
    Gene4 = B and C

More complex pattern are applied to our experiments to achieve detection of up to 60 single molecules....
"""