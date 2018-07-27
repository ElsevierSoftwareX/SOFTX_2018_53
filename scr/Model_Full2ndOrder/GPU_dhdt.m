function result = GPU_dhdt(dataqx,datapz)

result=-real(dataqx)-real(datapz)-imag(dataqx)*1i-imag(datapz)*1i;
end