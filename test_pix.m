function inside = test_pix(pix, absc_min, absc_max, ord_min, ord_max)

if(pix(1) < absc_min)
	inside = -1;
	return;
end
if(pix(2) < ord_min)
	inside = -1;
	return;
end
if(pix(1) > absc_max)
	inside = -1;
	return;
end
if(pix(2) > ord_max)
	inside = -1;
	return;
end
	
inside = 1;	
return;
