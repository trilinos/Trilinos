BEGIN {last = $2 }

/inuse/ { 
    
    if (last != $2)
    {
	print "\n"
    }
    print $0
    last = $2
}

