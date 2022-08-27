#!/usr/bin/python3

def subsetMerList(string, length):
    # Define this recursively to avoid repetition:
    if length == 1:
        # The base case is: Return the supplied string with A,C,G,T
        # tacked on the end, one element in the list for each.
        return [string + "A", string + "C", string + "G", string + "T"]
    else:
        # The recursive case is: Tack one of the letters onto the string
        # we got, and then call ourselves to add another layer, adding
        # the lists together.
        return subsetMerList(string + "A", length-1) + \
            subsetMerList(string + "C", length-1) + \
            subsetMerList(string + "G", length-1) + \
            subsetMerList(string + "T", length-1)


if __name__ == "__main__":
    for mer in subsetMerList("",4):
        print(mer + ".mers")


