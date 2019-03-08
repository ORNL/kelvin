Adding Constitutive Relationships
=

Constitutive Relationships/Equations can be added to Kelvin through a simple process:
1) Provide a class that implements the ConstitutiveRelationship interface found in ConstitutiveRelationship.h.
2) Register the new Constitutive Relationship with the ConstitutiveRelationshipService class in the configureConstitutiveRelationships() function in Kelvin.cpp.
3) Assign the appropriate material id in the particles file for all affected particles.

Implementing ConstitutiveRelationship
==

The ConstitutiveRelationship interface is defined in ConstitutiveRelationship.h. This class can be implemented by creating a subclass of an existing implementation, or starting from scratch. For a new constitutive relationship - called AwesomeRelationship for simplicity - create separate AwesomeRelationship.h and Awesome.cpp relationships in the source directory, as well as AwesomeRelationshipTest.cpp in the test directory. The src/HydrostaticCR* files are good examples of how to do this.

Registering with the Constitutive Relationship Service
==

Constitutive relationships are registered in the configureConstitutiveRelationships() function in Kelvin.cpp. The following code snippet is an example of how this registration is performed:
```cpp
	// Create the hydrostatic constitutive relationships
	unique_ptr<ConstitutiveRelationship> hydrostaticCR =
			make_unique<HydrostaticCR>(data);
	int hydrostaticID = 2;
	// Add it to the service
	ConstitutiveRelationshipService::add(hydrostaticID,std::move(hydrostaticCR));
```
This code works by creating a pointer to the constitutive relationship and then registering that pointer against a key - the material id - in the ConstitutiveRelationshipService. This service is a map that stores all of the available constitutive relationships so that they can be retrieved by their material id later.

Setting particle material ids
==

The particle material id is an integer that uniquely identifies a specific constitutive relationship and thus material type. In principle this id could be set by looping over all the particles in the code with either basic or sophisticated logic, but the most practical way to set the value is by editing the particles.csv file that describes the particles that will be used in the system. For example, the following would describe a particle with a position and material in particles.csv:
```
2.56728, 1.8767, 1.19445, 1
```
Modifying this line by changing the trailing "1" to a trailing "2" changes the material type for that particle accordingly from 1 to 2, as shown below:
```
2.56728, 1.8767, 1.19445, 2
```