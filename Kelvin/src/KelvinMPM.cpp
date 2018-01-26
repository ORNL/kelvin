
// Gets initial MPs from mesh

	auto & feSpace = meshContainer.getSpace();
	int numElements = feSpace.GetNE();
	for (int i = 0; i < numElements; i++) {
		auto * element = feSpace.GetFE(i);
		auto * transform = feSpace.GetElementTransformation(i);
		Vector point;
		auto & intRule = element->GetNodes();
		int numIntPoints = intRule.GetNPoints();
		for (int j = 0; j < numIntPoints; j++) {
			auto & intPoint = intRule.IntPoint(j);
			transform->Transform(intPoint,point);
			cout << point(0) << " " << point(1) << " " << point(2) << endl;
		}
	}
